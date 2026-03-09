"""
Responsibilities :
  - Loading the reference star catalog from the DB (read-only, no gap-filling)
  - Building astropy SkyCoords from catalog rows using J2000 coordinates
  - Computing Roman's observable windows via keepout
  - Checking solar angle and pitch angle constraints day-by-day per window
  - Returning ALL valid reference stars per window, ranked by observable days

OVERALL FLOW
1. User specifies band (1, 3, or 4) AND contrast level ('high' or 'med').
   These together select the correct grade column from the DB:
       Band 1  high  ==  st_psfgrade_nfb1_high
       Band 1  med   ==  st_psfgrade_nfb1_med
       Band 3  high  ==  st_psfgrade_specb3_high
       Band 3  med   ==  st_psfgrade_specb3_med
       Band 4  high  ==  st_psfgrade_wfb4_high
       Band 4  med   ==  st_psfgrade_wfb4_med

2. compute_keepout() on the science target.

3. For each observable window, check all candidate reference stars ranked by:
       a. Grade (A is best, then B, then C)
       b. Closest magnitude to the science target in the chosen band (V or I).
          If no science target magnitude is available, sort brightest-first.

4. For each candidate, check TWO constraints DAY BY DAY:
       a. Solar angle: can Roman physically point at the reference star
          (solar angle in [54°, 126°]) on that day?
       b. Pitch angle: is the pitch angle difference between science target
          and reference star less than 5 degrees on that day?

   A star is valid if it passes BOTH constraints on at least one day.
   All valid stars are returned, sorted by number of valid days (descending).
"""

import os
import numpy as np
import astropy.units as u
from astropy.time import Time
import astropy.coordinates as c
import pandas as pd
import sqlalchemy as sql

from roman_pointing.roman_observability import (
    get_target_coords,
    compute_roman_angles,
    compute_keepout,
)

GRADE_COLUMNS = {
    (1, 'high'): 'st_psfgrade_nfb1_high',
    (1, 'med'):  'st_psfgrade_nfb1_med',
    (3, 'high'): 'st_psfgrade_specb3_high',
    (3, 'med'):  'st_psfgrade_specb3_med',
    (4, 'high'): 'st_psfgrade_wfb4_high',
    (4, 'med'):  'st_psfgrade_wfb4_med',
}
ALL_GRADE_COLUMNS = list(GRADE_COLUMNS.values())

REF_GRADES = ['A', 'B', 'C']
SKIP_NAMES = {'-', 'TBD', '?', ''}

SUN_MIN        = 54    # minimum solar angle (deg)
SUN_MAX        = 126   # maximum solar angle (deg)
MAX_PITCH_DIFF = 5.0   # maximum pitch angle difference sci <-> ref (deg)

# Columns read from the Stars table
_CATALOG_COLUMNS = [
    'main_id', 'st_name',
    'ra', 'dec',
    'sy_vmag', 'sy_imag',
    'sy_dist', 'sy_plx',
    'sy_pmra', 'sy_pmdec',
    'st_radv', 'spectype',
]

def load_catalog(engine: sql.engine.base.Engine) -> pd.DataFrame:
    """Load the reference star catalog from the DB. Read-only; no gap-filling.

    Only stars that carry a grade (A/B/C) in at least one of the six grade
    columns are returned.  Numeric columns are coerced; nothing is written
    back to the database.

    Args:
        engine: SQLAlchemy engine connected to plandb.

    Returns:
        pd.DataFrame with one row per graded star and all columns needed
        by select_ref_star().

    Raises:
        RuntimeError: if no populated grade columns are found.
    """
    metadata    = sql.MetaData()
    stars_table = sql.Table('Stars', metadata, autoload_with=engine)
    actual_cols = {col.name for col in stars_table.columns}

    # Determine which grade columns exist and are populated
    present_grade_cols = [col for col in ALL_GRADE_COLUMNS if col in actual_cols]

    all_null = True
    if present_grade_cols:
        with engine.connect() as conn:
            for gcol in present_grade_cols:
                count = conn.execute(
                    sql.text(f"SELECT COUNT(*) FROM Stars WHERE {gcol} IS NOT NULL")
                ).scalar()
                if count and count > 0:
                    all_null = False
                    break

    if not present_grade_cols or all_null:
        fallback = 'st_psfgrade'
        if fallback in actual_cols:
            print(f"  New grade columns are all NULL — falling back to '{fallback}'.")
            present_grade_cols = [fallback]
        else:
            raise RuntimeError(
                f"No populated grade columns found in Stars table. "
                f"Tried {ALL_GRADE_COLUMNS} (all NULL) and '{fallback}' (absent). "
                f"Available columns: {sorted(actual_cols)}"
            )

    absent = set(ALL_GRADE_COLUMNS) - set(present_grade_cols)
    if absent:
        print(f"  Note: grade column(s) not in DB yet: {sorted(absent)}")

    # Select only graded stars
    grade_sel = [
        getattr(stars_table.c, col)
        for col in present_grade_cols
        if col in actual_cols
    ]
    base_sel = [
        getattr(stars_table.c, col)
        for col in _CATALOG_COLUMNS
        if col in actual_cols
    ]

    # WHERE: at least one grade column is non-NULL
    where_clause = sql.or_(
        *[getattr(stars_table.c, col).isnot(None) for col in present_grade_cols]
    )

    stmt = sql.select(*base_sel, *grade_sel).where(where_clause)

    with engine.connect() as conn:
        df = pd.read_sql(stmt, conn)

    # Add placeholder columns for grade columns absent from this DB version
    for col in ALL_GRADE_COLUMNS:
        if col not in df.columns:
            df[col] = np.nan

    # Coerce numeric columns
    for col in ('sy_vmag', 'sy_imag', 'sy_dist', 'sy_plx',
                'sy_pmra', 'sy_pmdec', 'st_radv'):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Derive sy_dist from parallax in-memory (not written to DB)
    missing_dist = df['sy_dist'].isna() & df['sy_plx'].notna() & (df['sy_plx'] > 0)
    n_derived = int(missing_dist.sum())
    if n_derived:
        from astropy.coordinates import Distance
        df.loc[missing_dist, 'sy_dist'] = Distance(
            parallax=df.loc[missing_dist, 'sy_plx'].values * u.mas
        ).pc
        print(f"  Derived sy_dist from sy_plx for {n_derived} star(s) (in memory only).")

    # Convenience columns used by select_ref_star
    df['mag_v'] = df['sy_vmag']
    df['mag_i'] = df['sy_imag']

    print(f"Catalog loaded: {len(df)} graded reference star(s).")
    return df


def get_science_mag(sci_name: str, band: int, engine=None) -> float | None:
    """Look up the science target magnitude in the Stars table (read-only).

    Tries main_id and st_name.  Returns None if not found or if the value
    is NULL — the caller gracefully falls back to brightest-first sorting.

    Args:
        sci_name: Science target name (SIMBAD-resolvable).
        band:     1 == V band (sy_vmag), 3 or 4 ==I band (sy_imag).
        engine:   SQLAlchemy engine.  If None, returns None immediately.

    Returns:
        float magnitude, or None.
    """
    if engine is None:
        return None

    mag_col  = 'sy_vmag' if band == 1 else 'sy_imag'
    band_lbl = 'V' if band == 1 else 'I'

    try:
        metadata    = sql.MetaData()
        stars_table = sql.Table('Stars', metadata, autoload_with=engine)
        stmt = sql.select(getattr(stars_table.c, mag_col)).where(
            sql.or_(
                stars_table.c.main_id == sci_name,
                stars_table.c.st_name == sci_name,
            )
        ).limit(1)

        with engine.connect() as conn:
            row = conn.execute(stmt).fetchone()

        if row is not None and row[0] is not None:
            val = float(row[0])
            if not np.isnan(val):
                print(f"  Science target {band_lbl}-band mag: {val:.2f}")
                return val
    except Exception as e:
        print(f"  Warning: could not look up science magnitude ({e}).")

    print(f"  Science target {band_lbl}-band mag not found in DB; "
          f"will sort brightest-first.")
    return None

def build_skycoord(star) -> c.SkyCoord:
    """Build an astropy SkyCoord in BarycentricMeanEcliptic from a star record.

    Uses J2000 coordinates (ra, dec). Proper motion, parallax, and radial
    velocity are included when available.

    Args:
        star: dict or pd.Series with astrometric columns:
              ra, dec, sy_plx, sy_dist, sy_pmra, sy_pmdec, st_radv.

    Returns:
        SkyCoord in BarycentricMeanEcliptic.
    """
    def val(key, fallback=None):
        v = star[key] if isinstance(star, dict) else star.get(key, fallback)
        return (
            None if (v is None or (isinstance(v, float) and np.isnan(v)))
            else float(v)
        )

    kwargs = dict(
        ra=val('ra') * u.degree,
        dec=val('dec') * u.degree,
        frame='icrs',
        equinox='J2000',
        obstime='J2000',
    )

    if val('sy_plx'):    kwargs['distance']        = c.Distance(parallax=val('sy_plx') * u.mas)
    elif val('sy_dist'): kwargs['distance']        = val('sy_dist') * u.parsec
    if val('sy_pmra'):   kwargs['pm_ra_cosdec']    = val('sy_pmra')  * u.mas / u.yr
    if val('sy_pmdec'):  kwargs['pm_dec']          = val('sy_pmdec') * u.mas / u.yr
    if val('st_radv'):   kwargs['radial_velocity'] = val('st_radv')  * u.km / u.s

    return c.SkyCoord(**kwargs).transform_to(c.BarycentricMeanEcliptic)


def get_observable_windows(ts: Time, keepout_array: np.ndarray) -> list:
    """Extract contiguous time windows when Roman can observe the science target.

    Args:
        ts (astropy.time.Time):     Array of time values matching keepout_array.
        keepout_array (np.ndarray): Boolean array; True = observable.

    Returns:
        list of tuples:
            (start_Time, end_Time, start_str, end_str, duration_days)
    """
    windows   = []
    in_window = False
    start_idx = 0

    for i, obs in enumerate(keepout_array):
        if obs and not in_window:
            in_window = True
            start_idx = i
        elif not obs and in_window:
            in_window = False
            windows.append((
                ts[start_idx], ts[i - 1],
                ts[start_idx].iso.split('T')[0],
                ts[i - 1].iso.split('T')[0],
                ts[i - 1].mjd - ts[start_idx].mjd,
            ))

    if in_window:
        windows.append((
            ts[start_idx], ts[-1],
            ts[start_idx].iso.split('T')[0],
            ts[-1].iso.split('T')[0],
            ts[-1].mjd - ts[start_idx].mjd,
        ))

    return windows


def check_ref_in_window(
    ref_coord: c.SkyCoord,
    ref_name: str,
    win_start: Time,
    win_end: Time,
    sci_pitch_in_window: np.ndarray,
) -> tuple:
    """Check how many days a reference star satisfies both pointing constraints.

    Checks constraints day-by-day. A day passes if both the solar angle
    and pitch angle constraints are met on that specific day. The star is
    considered valid if it passes on at least one day within the window.

    Args:
        ref_coord:           SkyCoord of the reference star (BarycentricMeanEcliptic)
        ref_name:            Name string (for logging)
        win_start:           Start of the observable window
        win_end:             End of the observable window
        sci_pitch_in_window: Science target pitch angles over this window (deg)

    Returns:
        tuple[bool, int, float, np.ndarray]:
            passes            -- True if valid on at least one day
            n_days_valid      -- Number of days passing both constraints
            min_pitch_diff    -- Minimum pitch difference on any valid day (deg);
                                 999 if no valid days
            pitch_diff_series -- Per-day |delta_pitch| (deg), NaN where solar
                                 angle blocks observation
    """
    duration_days = win_end.mjd - win_start.mjd
    if duration_days <= 0:
        return False, 0, 999.0, np.array([])

    start_str = win_start.isot if hasattr(win_start, 'isot') else str(win_start)

    _, ref_sun_ang, _, ref_pitch = compute_roman_angles(
        ref_coord, start_str, duration_days, time_step=1.0
    )
    ref_sun_d   = ref_sun_ang.to(u.degree).value
    ref_pitch_d = ref_pitch.to(u.degree).value

    min_len    = min(len(sci_pitch_in_window), len(ref_sun_d), len(ref_pitch_d))
    solar_ok   = (ref_sun_d[:min_len] > SUN_MIN) & (ref_sun_d[:min_len] < SUN_MAX)
    pitch_diff = np.abs(sci_pitch_in_window[:min_len] - ref_pitch_d[:min_len])
    pitch_ok   = pitch_diff < MAX_PITCH_DIFF

    pitch_diff_series            = pitch_diff.copy().astype(float)
    pitch_diff_series[~solar_ok] = np.nan   # NaN = solar angle blocked

    valid_days = solar_ok & pitch_ok
    n_valid    = int(np.sum(valid_days))

    if n_valid == 0:
        return False, 0, 999.0, pitch_diff_series

    return True, n_valid, float(np.min(pitch_diff[valid_days])), pitch_diff_series

def select_ref_star(
    sci_name: str,
    analysis_start: str,
    analysis_days: float,
    band: int,
    contrast: str,
    catalog: pd.DataFrame,
    engine=None,
    time_step: float = 1.0,
) -> dict:
    """Find all valid reference stars for each observable window of a science target.

    Reference stars are checked day-by-day within each window. A star is
    valid if it passes both the solar angle and pitch angle constraints on
    at least one day. Results are sorted by number of valid days (descending).

    Args:
        sci_name (str):         Science target name recognized by SIMBAD.
        analysis_start (str):   ISO start date, e.g. '2027-01-01T00:00:00'.
        analysis_days (float):  Total span to analyze (days).
        band (int):             1 = V band (NFB), 3 = spec I band, 4 = wide I band.
        contrast (str):         'high' or 'med'.
        catalog (pd.DataFrame): Reference star catalog from load_catalog().
        engine:                 SQLAlchemy engine for plandb (optional;
                                used only for science target magnitude lookup).
        time_step (float):      Time resolution for angle calculations (days).

    Returns:
        dict with keys:
            'science_target', 'band', 'contrast', 'grade_column',
            'visibility_pct', 'sort_method',
            'observable_windows' (list of per-window result dicts), each with:
                'start', 'end', 'duration_days',
                'valid_refs' (list of all passing ref stars, sorted by n_valid_days),
                'best_ref'   (ref star with most valid days, or None),
                'pitch_df'   (DataFrame: rows=dates, cols=valid ref star names),
                'pitch_csv'  (path to saved CSV).

    Raises:
        ValueError: if (band, contrast) is not a recognised combination.
        ValueError: if neither a pre-resolved 'grade' column nor the raw grade
                    column is present in the catalog.
    """
    key = (band, contrast.lower())
    if key not in GRADE_COLUMNS:
        valid = ', '.join(f"band={b} contrast={ct}" for b, ct in GRADE_COLUMNS)
        raise ValueError(
            f"Unknown (band={band}, contrast='{contrast}'). "
            f"Valid combinations: {valid}"
        )

    grade_col  = GRADE_COLUMNS[key]
    mag_col    = 'mag_v' if band == 1 else 'mag_i'
    band_lbl   = 'V'     if band == 1 else 'I'

    print(f"\nUsing grade column: {grade_col}  |  mag column: {mag_col}")

    # Resolve grade into a plain 'grade' column.
    # catalog.py may have already done this (column is named 'grade').
    # If not, fall back to reading the raw DB column name.
    candidates = catalog.copy()
    if 'grade' in candidates.columns:
        candidates['grade'] = candidates['grade'].str.strip()
    elif grade_col in candidates.columns:
        candidates['grade'] = candidates[grade_col].str.strip()
    else:
        raise ValueError(
            f"Neither a pre-resolved 'grade' column nor the raw grade column "
            f"'{grade_col}' was found in the catalog. "
            f"Available columns: {list(catalog.columns)}"
        )
    candidates = candidates[candidates['grade'].isin(REF_GRADES)].copy()

    # Science target coords
    print(f"Querying SIMBAD for science target '{sci_name}'...")
    coords = get_target_coords([sci_name])
    if sci_name not in coords:
        return {'error': f"Science target '{sci_name}' not found in SIMBAD."}
    sci_coord = coords[sci_name]
    print(f"  Found '{sci_name}'.\n")

    # Science target magnitude (DB lookup only — no external queries)
    print(f"Looking up {band_lbl}-band magnitude for '{sci_name}'...")
    sci_mag = get_science_mag(sci_name, band, engine=engine)

    # Keepout / observable windows
    print(f"\nComputing Roman visibility for '{sci_name}' "
          f"over {analysis_days:.0f} days...")
    ts, keepout, _ = compute_keepout(
        {sci_name: sci_coord}, analysis_start, analysis_days, time_step
    )
    sci_keepout    = keepout[sci_name]
    visibility_pct = (np.sum(sci_keepout) / len(sci_keepout)) * 100
    print(f"  Observable {visibility_pct:.1f}% of the time.")

    windows = get_observable_windows(ts, sci_keepout)
    if not windows:
        return {
            'science_target': sci_name,
            'band':           band,
            'contrast':       contrast,
            'grade_column':   grade_col,
            'error': (f"'{sci_name}' is never observable by Roman "
                      f"during this period."),
        }

    print(f"\nFound {len(windows)} observable window(s):")
    for i, (_, _, ws, we, wd) in enumerate(windows):
        print(f"  Window {i + 1}: {ws} to {we} ({wd:.1f} days)")

    # Science pitch angles (computed once, sliced per window)
    _, _, _, sci_pitch_full = compute_roman_angles(
        sci_coord, analysis_start, analysis_days, time_step
    )
    sci_pitch_vals = sci_pitch_full.to(u.degree).value

    # Rank candidates by grade then magnitude
    candidates = candidates.dropna(subset=[mag_col])
    grade_order = {g: i for i, g in enumerate(REF_GRADES)}
    candidates['grade_rank'] = candidates['grade'].map(grade_order).fillna(99)

    if sci_mag is not None:
        candidates['mag_diff'] = (candidates[mag_col] - sci_mag).abs()
        candidates  = candidates.sort_values(['grade_rank', 'mag_diff'])
        sort_method = (f"grade (A>B>C) then closest {band_lbl}-band "
                       f"magnitude to science target")
    else:
        candidates['mag_diff'] = np.nan
        candidates  = candidates.sort_values(['grade_rank', mag_col])
        sort_method = f"grade (A>B>C) then brightest {band_lbl}"

    print(f"\nSorting by: {sort_method}")
    print(f"Candidates available: {len(candidates)}")

    # Build SkyCoords (J2000 only)
    print("Building reference star coordinates (J2000)...")
    ref_coords = {}
    for _, ref in candidates.iterrows():
        name = ref['main_id']
        if not isinstance(name, str) or name.strip() in SKIP_NAMES:
            continue
        ref_coords[name] = build_skycoord(ref)

    print(f"  Successfully built coordinates for {len(ref_coords)} stars.")

    # Search each window
    results = []
    print("\nSearching for reference stars in each window...")

    for win_idx, (win_start, win_end, ws, we, wd) in enumerate(windows):
        print(f"\nWindow {win_idx + 1}: {ws} to {we} ({wd:.1f} days)")

        win_start_idx = int((win_start.mjd - ts[0].mjd) / time_step)
        win_end_idx   = int((win_end.mjd   - ts[0].mjd) / time_step)
        sci_pitch_win = sci_pitch_vals[win_start_idx:win_end_idx + 1]

        valid_refs   = []
        pitch_series = {}   # ref_name -> per-day |Δpitch| array

        for grade in REF_GRADES:
            grade_cands = candidates[candidates['grade'] == grade]
            if grade_cands.empty:
                continue

            for _, ref in grade_cands.iterrows():
                ref_name = ref['main_id']
                if ref_name not in ref_coords:
                    continue

                passes, n_days, min_pitch, pd_series = check_ref_in_window(
                    ref_coords[ref_name], ref_name,
                    win_start, win_end, sci_pitch_win,
                )
                pitch_series[ref_name] = pd_series

                if passes:
                    valid_refs.append({
                        'reference_star': ref_name,
                        'grade':          grade,
                        'mag_diff':       None if np.isnan(ref['mag_diff']) else float(ref['mag_diff']),
                        'n_valid_days':   n_days,
                        'min_pitch_diff': min_pitch,
                    })

        # Sort by most valid days first
        valid_refs.sort(key=lambda r: r['n_valid_days'], reverse=True)

        # Build pitch DataFrame: rows=dates, cols=ref stars (valid only)
        n_days_win = int(round(win_end.mjd - win_start.mjd)) + 1
        dates = [
            (win_start + i * u.day).to_value('iso', subfmt='date')
            for i in range(n_days_win)
        ]
        valid_names = [r['reference_star'] for r in valid_refs]
        pitch_data  = {}
        for name in valid_names:
            series = pitch_series.get(name, np.array([]))
            vals   = list(series[:len(dates)])
            vals  += [np.nan] * (len(dates) - len(vals))
            pitch_data[name] = vals

        pitch_df            = pd.DataFrame(pitch_data, index=dates)
        pitch_df.index.name = 'date'

        # Save pitch table CSV
        safe_name = sci_name.replace(' ', '_').replace('*', '').strip('_')
        csv_path  = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            f"pitch_table_{safe_name}_band{band}_{contrast}_window{win_idx + 1}.csv",
        )
        pitch_df.to_csv(csv_path, float_format='%.3f')
        print(f"  Pitch table saved → {csv_path}")

        print(f"  Found {len(valid_refs)} valid reference star(s) for this window.")
        results.append({
            'start':         ws,
            'end':           we,
            'duration_days': wd,
            'valid_refs':    valid_refs,
            'best_ref':      valid_refs[0] if valid_refs else None,
            'pitch_df':      pitch_df,
            'pitch_csv':     csv_path,
        })

    return {
        'science_target':     sci_name,
        'band':               band,
        'contrast':           contrast,
        'grade_column':       grade_col,
        'visibility_pct':     visibility_pct,
        'sort_method':        sort_method,
        'observable_windows': results,
    }

if __name__ == "__main__":
    import corgidb.ingest

    print("ReferenceStarPicker\n")

    SCIENCE_TARGET = "47 Uma"
    BAND           = 1        # 1 = NFB (V-band), 3 = spec (I-band), 4 = wide (I-band)
    CONTRAST       = 'high'   # 'high' or 'med' — combined with BAND selects grade column:
                              #   band=1 high == st_psfgrade_nfb1_high
                              #   band=1 med  == st_psfgrade_nfb1_med
                              #   band=3 high == st_psfgrade_specb3_high
                              #   band=3 med  == st_psfgrade_specb3_med
                              #   band=4 high == st_psfgrade_wfb4_high
                              #   band=4 med  == st_psfgrade_wfb4_med
    ANALYSIS_START = "2027-01-01T00:00:00"
    ANALYSIS_DAYS  = 365
 

    eng = corgidb.ingest.gen_engine('plandb_user', 'plandb_scratch')

    catalog = load_catalog(eng)
    print(f"Missing mag_v: {catalog['mag_v'].isna().sum()}")
    print(f"Missing mag_i: {catalog['mag_i'].isna().sum()}")
    print(f"Grade value counts:\n{catalog['st_psfgrade_nfb1_high'].value_counts(dropna=False)}")
    print(f"\nCatalog ready: {len(catalog)} reference stars.\n")

    result = select_ref_star(
        SCIENCE_TARGET, ANALYSIS_START, ANALYSIS_DAYS,
        band=BAND, contrast=CONTRAST,
        catalog=catalog, engine=eng,
    )

    print("\n" + "=" * 60)
    print(f"RESULTS FOR: {result.get('science_target')} "
          f"(Band {result.get('band')}, {result.get('contrast')} contrast)")
    print(f"Grade column: {result.get('grade_column')}")

    if 'error' in result:
        print(f"ERROR: {result['error']}")
    else:
        print(f"Observable {result['visibility_pct']:.1f}% of the time")
        print(f"Sort method: {result['sort_method']}")
        print("=" * 60)
        for i, win in enumerate(result['observable_windows']):
            print(f"\nWindow {i + 1}: {win['start']} to {win['end']} "
                  f"({win['duration_days']:.1f} days)")
            if not win['valid_refs']:
                print("  No suitable reference stars found.")
            else:
                print(f"  {len(win['valid_refs'])} valid reference star(s) "
                      f"(sorted by most observable days):")
                for ref in win['valid_refs']:
                    mag_str = (
                        f", mag diff={ref['mag_diff']:.2f}"
                        if ref['mag_diff'] is not None else ""
                    )
                    print(f"    {ref['reference_star']} "
                          f"(grade={ref['grade']}{mag_str}, "
                          f"{ref['n_valid_days']} valid day(s), "
                          f"min pitch diff={ref['min_pitch_diff']:.4f} deg)")