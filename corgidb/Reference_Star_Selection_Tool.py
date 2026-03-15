"""Reference star selection module for the Roman Space Telescope.

    - Loading the reference star catalog from the DB (read-only, no gap-filling).
    - Building astropy SkyCoords from catalog rows using J2000 coordinates.
    - Computing Roman's observable windows via keepout.
    - Checking solar angle and pitch angle constraints day-by-day per window.
    - Returning all valid reference stars per window, ranked by grade then metric.

Overall Flow:
    1. User specifies band (1, 3, or 4) and contrast level ('high' or 'med').
       These together select the correct grade column from the DB::

           Band 1  high  ==  st_psfgrade_nfb1_high
           Band 1  med   ==  st_psfgrade_nfb1_med
           Band 3  high  ==  st_psfgrade_specb3_high
           Band 3  med   ==  st_psfgrade_specb3_med
           Band 4  high  ==  st_psfgrade_wfb4_high
           Band 4  med   ==  st_psfgrade_wfb4_med

    2. ``compute_keepout()`` is called on the science target.

    3. For each observable window, all passing reference stars are sorted by:

           PRIMARY   — grade tier (order from allowed_grades; A before B before C
                       by default). This is always the outermost ranking level.

           SECONDARY — the chosen sort_mode, applied purely within each grade tier:

               'valid_days'      most valid days first (descending n_valid_days)
               'closest_mag'     closest |ref_mag - sci_mag| first (ascending);
                                  falls back to 'brightest' when sci_mag unavailable
               'brightest'       brightest first (ascending magnitude)
               'faintest'        faintest first (descending magnitude)
               'closest_pitch'   smallest min_pitch_diff first (ascending)
               'farthest_pitch'  largest min_pitch_diff first (descending)

           There is no tertiary tiebreaker. Within each grade tier the secondary
           metric is the sole ranking criterion. Stars with identical metric
           values retain their iteration order (stable sort).

    4. For each candidate, two constraints are checked day by day:

           a. Solar angle: Roman can point at the reference star
              (solar angle in [54 deg, 126 deg]) on that day.
           b. Pitch angle: pitch difference between science target and reference
              star is less than 5 degrees on that day.

       A star is valid if it passes both constraints on at least one day.
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

SUN_MIN = 54
SUN_MAX = 126
MAX_PITCH_DIFF = 5.0

SORT_MODES = (
    'valid_days',
    'closest_mag',
    'brightest',
    'faintest',
    'closest_pitch',
    'farthest_pitch',
)

SORT_MODE_LABELS = {
    'valid_days':     'Most valid days first',
    'closest_mag':    'Closest magnitude to science target',
    'brightest':      'Brightest first (ascending magnitude)',
    'faintest':       'Faintest first (descending magnitude)',
    'closest_pitch':  'Smallest pitch angle difference first',
    'farthest_pitch': 'Largest pitch angle difference first',
}

CATALOG_COLUMNS = [
    'main_id', 'st_name',
    'ra', 'dec',
    'sy_vmag', 'sy_imag',
    'sy_dist', 'sy_plx',
    'sy_pmra', 'sy_pmdec',
    'st_radv', 'spectype',
]

LARGE_SENTINEL = 1e9


def safe_float(value):
    """Return the float representation of a value, or None if missing or NaN.

    Args:
        value: Any scalar value to convert.

    Returns:
        float or None: The converted float, or None if conversion fails or
        the value is NaN.
    """
    if value is None:
        return None
    try:
        result = float(value)
        return None if np.isnan(result) else result
    except (TypeError, ValueError):
        return None


def make_sort_key(sort_mode):
    """Return a key function for sorting a list of valid-reference-star dicts.

    Sorting contract:
        PRIMARY  : ``grade_rank`` (int, lower = better grade, always first).
        SECONDARY: the chosen metric, computed purely from the ref dict fields.
                   No additional tiebreakers. The sort is stable so equal-metric
                   stars retain their original iteration order.

    Each returned key is a 2-tuple: ``(grade_rank, metric_value)``.
    An ascending sort on this tuple gives the correct final order.

    Args:
        sort_mode (str): One of the values in ``SORT_MODES``.

    Returns:
        callable: A key function suitable for use with ``list.sort(key=...)``.

    Raises:
        ValueError: If ``sort_mode`` is not a recognised value.
    """
    if sort_mode == 'valid_days':
        def key(ref):
            return (ref['grade_rank'], -ref['n_valid_days'])

    elif sort_mode == 'closest_mag':
        def key(ref):
            diff = ref['mag_diff'] if ref['mag_diff'] is not None else LARGE_SENTINEL
            return (ref['grade_rank'], diff)

    elif sort_mode == 'brightest':
        def key(ref):
            mag = ref['mag'] if ref['mag'] is not None else LARGE_SENTINEL
            return (ref['grade_rank'], mag)

    elif sort_mode == 'faintest':
        def key(ref):
            mag = ref['mag'] if ref['mag'] is not None else -LARGE_SENTINEL
            return (ref['grade_rank'], -mag)

    elif sort_mode == 'closest_pitch':
        def key(ref):
            pitch = ref['min_pitch_diff'] if ref['min_pitch_diff'] < 999 else LARGE_SENTINEL
            return (ref['grade_rank'], pitch)

    elif sort_mode == 'farthest_pitch':
        def key(ref):
            pitch = ref['min_pitch_diff'] if ref['min_pitch_diff'] < 999 else -LARGE_SENTINEL
            return (ref['grade_rank'], -pitch)

    else:
        raise ValueError(
            f"Unknown sort_mode='{sort_mode}'. Valid options: {SORT_MODES}"
        )

    return key


def load_catalog(engine):
    """Load the reference star catalog from the database.

    Reads are performed in read-only fashion; no gap-filling or writes occur.
    Only stars that carry a grade (A, B, or C) in at least one of the six
    grade columns are returned. Numeric columns are coerced in memory.

    Args:
        engine (sqlalchemy.engine.base.Engine): SQLAlchemy engine connected
            to plandb.

    Returns:
        pandas.DataFrame: One row per graded star.

    Raises:
        RuntimeError: If no populated grade columns are found in the database.
    """
    metadata = sql.MetaData()
    stars_table = sql.Table('Stars', metadata, autoload_with=engine)
    actual_cols = {col.name for col in stars_table.columns}

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

    grade_sel = [
        getattr(stars_table.c, col)
        for col in present_grade_cols
        if col in actual_cols
    ]
    base_sel = [
        getattr(stars_table.c, col)
        for col in CATALOG_COLUMNS
        if col in actual_cols
    ]

    where_clause = sql.or_(
        *[getattr(stars_table.c, col).isnot(None) for col in present_grade_cols]
    )
    stmt = sql.select(*base_sel, *grade_sel).where(where_clause)

    with engine.connect() as conn:
        df = pd.read_sql(stmt, conn)

    for col in ALL_GRADE_COLUMNS:
        if col not in df.columns:
            df[col] = np.nan

    for col in ('sy_vmag', 'sy_imag', 'sy_dist', 'sy_plx',
                'sy_pmra', 'sy_pmdec', 'st_radv'):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    missing_dist = (
        df['sy_dist'].isna()
        & df['sy_plx'].notna()
        & (df['sy_plx'] > 0)
    )
    n_derived = int(missing_dist.sum())
    if n_derived:
        from astropy.coordinates import Distance
        df.loc[missing_dist, 'sy_dist'] = Distance(
            parallax=df.loc[missing_dist, 'sy_plx'].values * u.mas
        ).pc
        print(
            f"  Derived sy_dist from sy_plx for {n_derived} star(s) (in memory only)."
        )

    df['mag_v'] = df['sy_vmag']
    df['mag_i'] = df['sy_imag']

    print(f"Catalog loaded: {len(df)} graded reference star(s).")
    return df


def get_science_mag(sci_name, band, engine=None):
    """Look up the science target magnitude in the Stars table.

    The database is accessed in read-only fashion. Returns ``None`` if the
    magnitude is not found, in which case the caller falls back to
    brightest-first sorting.

    Args:
        sci_name (str): SIMBAD-resolvable name of the science target.
        band (int): Photometric band. 1 selects V-band; any other value
            selects I-band.
        engine (sqlalchemy.engine.base.Engine, optional): SQLAlchemy engine.
            Defaults to None, which causes the function to return None
            immediately.

    Returns:
        float or None: The magnitude value, or None if not found.
    """
    if engine is None:
        return None

    mag_col = 'sy_vmag' if band == 1 else 'sy_imag'
    band_label = 'V' if band == 1 else 'I'

    try:
        metadata = sql.MetaData()
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
                print(f"  Science target {band_label}-band mag: {val:.2f}")
                return val
    except Exception as exc:
        print(f"  Warning: could not look up science magnitude ({exc}).")

    print(
        f"  Science target {band_label}-band mag not found in DB; "
        f"will sort brightest-first."
    )
    return None


def build_skycoord(star):
    """Build an astropy SkyCoord in BarycentricMeanEcliptic from a star record.

    Proper motion, parallax, distance, and radial velocity are included when
    available. The coordinate is first constructed in ICRS with J2000 equinox
    and then transformed to ``BarycentricMeanEcliptic``.

    Args:
        star (dict or pandas.Series): A mapping containing at minimum the keys
            ``'ra'`` and ``'dec'`` (in degrees, J2000). Optional keys are
            ``'sy_plx'`` (parallax in mas), ``'sy_dist'`` (distance in pc),
            ``'sy_pmra'`` (proper motion in RA, mas/yr),
            ``'sy_pmdec'`` (proper motion in Dec, mas/yr), and
            ``'st_radv'`` (radial velocity in km/s).

    Returns:
        astropy.coordinates.SkyCoord: The star's position in
        ``BarycentricMeanEcliptic``.
    """
    def val(key, fallback=None):
        raw = star[key] if isinstance(star, dict) else star.get(key, fallback)
        return (
            None
            if (raw is None or (isinstance(raw, float) and np.isnan(raw)))
            else float(raw)
        )

    kwargs = dict(
        ra=val('ra') * u.degree,
        dec=val('dec') * u.degree,
        frame='icrs',
        equinox='J2000',
        obstime='J2000',
    )

    if val('sy_plx'):
        kwargs['distance'] = c.Distance(parallax=val('sy_plx') * u.mas)
    elif val('sy_dist'):
        kwargs['distance'] = val('sy_dist') * u.parsec

    if val('sy_pmra'):
        kwargs['pm_ra_cosdec'] = val('sy_pmra') * u.mas / u.yr
    if val('sy_pmdec'):
        kwargs['pm_dec'] = val('sy_pmdec') * u.mas / u.yr
    if val('st_radv'):
        kwargs['radial_velocity'] = val('st_radv') * u.km / u.s

    return c.SkyCoord(**kwargs).transform_to(c.BarycentricMeanEcliptic)


def get_observable_windows(times, keepout_array):
    """Extract contiguous observable windows from a keepout boolean array.

    Args:
        times (astropy.time.Time): Array of times corresponding to each
            element of ``keepout_array``.
        keepout_array (numpy.ndarray): Boolean array where ``True`` indicates
            the target is observable.

    Returns:
        list of tuple: Each tuple contains:
            - start_time (astropy.time.Time)
            - end_time (astropy.time.Time)
            - start_str (str): ISO date string of the window start.
            - end_str (str): ISO date string of the window end.
            - duration_days (float): Length of the window in days.
    """
    windows = []
    in_window = False
    start_idx = 0

    for i, observable in enumerate(keepout_array):
        if observable and not in_window:
            in_window = True
            start_idx = i
        elif not observable and in_window:
            in_window = False
            windows.append((
                times[start_idx],
                times[i - 1],
                times[start_idx].iso.split('T')[0],
                times[i - 1].iso.split('T')[0],
                times[i - 1].mjd - times[start_idx].mjd,
            ))

    if in_window:
        windows.append((
            times[start_idx],
            times[-1],
            times[start_idx].iso.split('T')[0],
            times[-1].iso.split('T')[0],
            times[-1].mjd - times[start_idx].mjd,
        ))

    return windows


def check_ref_in_window(ref_coord, win_start, win_end, sci_pitch_in_window):
    """Check daily solar-angle and pitch-angle constraints for one reference star.

    For each day in the window the function evaluates two conditions:

        1. The reference star's solar angle is within ``[SUN_MIN, SUN_MAX]``
           degrees, meaning Roman can actually point at it.
        2. The absolute pitch-angle difference between the science target and
           the reference star is less than ``MAX_PITCH_DIFF`` degrees.

    Args:
        ref_coord (astropy.coordinates.SkyCoord): The reference star's
            position in ``BarycentricMeanEcliptic``.
        win_start (astropy.time.Time): Start of the observable window.
        win_end (astropy.time.Time): End of the observable window.
        sci_pitch_in_window (numpy.ndarray): Array of science-target pitch
            angles (degrees) for each day in the window.

    Returns:
        tuple: A 5-element tuple containing:
            - passes (bool): ``True`` if at least one day satisfies both
              constraints.
            - n_days_valid (int): Number of days satisfying both constraints.
            - min_pitch_diff (float): Minimum pitch difference across valid
              days; 999.0 if no valid days exist.
            - pitch_diff_series (numpy.ndarray): Pitch differences per day;
              days that fail the solar-angle check are set to NaN.
            - valid_mask (numpy.ndarray of bool): Per-day validity flags.
    """
    duration_days = win_end.mjd - win_start.mjd
    if duration_days <= 0:
        return False, 0, 999.0, np.array([]), np.array([], dtype=bool)

    start_str = win_start.isot if hasattr(win_start, 'isot') else str(win_start)

    _, ref_sun_ang, _, ref_pitch = compute_roman_angles(
        ref_coord, start_str, duration_days, time_step=1.0
    )
    ref_sun_d = ref_sun_ang.to(u.degree).value
    ref_pitch_d = ref_pitch.to(u.degree).value

    min_len = min(len(sci_pitch_in_window), len(ref_sun_d), len(ref_pitch_d))
    solar_ok = (ref_sun_d[:min_len] > SUN_MIN) & (ref_sun_d[:min_len] < SUN_MAX)
    pitch_diff = np.abs(sci_pitch_in_window[:min_len] - ref_pitch_d[:min_len])
    pitch_ok = pitch_diff < MAX_PITCH_DIFF

    pitch_diff_series = pitch_diff.copy().astype(float)
    pitch_diff_series[~solar_ok] = np.nan

    valid_mask = solar_ok & pitch_ok
    n_valid = int(np.sum(valid_mask))

    if n_valid == 0:
        return False, 0, 999.0, pitch_diff_series, valid_mask

    return (
        True,
        n_valid,
        float(np.min(pitch_diff[valid_mask])),
        pitch_diff_series,
        valid_mask,
    )


def select_ref_star(
    sci_name,
    analysis_start,
    analysis_days,
    band,
    contrast,
    catalog,
    engine=None,
    time_step=1.0,
    allowed_grades=None,
    sort_mode='valid_days',
):
    """Find all valid reference stars for each observable window of a science target.

    Grade tier is always the primary sort key. Within each grade tier the
    chosen ``sort_mode`` is applied as the sole secondary criterion:

        - ``valid_days``     descending ``n_valid_days``
        - ``closest_mag``     ascending  ``|ref_mag - sci_mag|``
        - ``brightest``       ascending  ``ref_mag``
        - ``faintest``        descending ``ref_mag``
        - ``closest_pitch``   ascending  ``min_pitch_diff``
        - ``farthest_pitch``  descending ``min_pitch_diff``

    Args:
        sci_name (str): SIMBAD-resolvable science target name.
        analysis_start (str): ISO start date, e.g. ``'2027-01-01T00:00:00'``.
        analysis_days (float): Total span to analyse in days.
        band (int): Photometric band. 1 = V-band (NFB), 3 = spec I-band,
            4 = wide I-band.
        contrast (str): ``'high'`` or ``'med'``.
        catalog (pandas.DataFrame): Reference star table from
            :func:`load_catalog`.
        engine (sqlalchemy.engine.base.Engine, optional): SQLAlchemy engine
            used for science-target magnitude lookup. Defaults to None.
        time_step (float, optional): Time resolution for angle calculations
            in days. Defaults to 1.0.
        allowed_grades (list of str, optional): Which grade tiers to include
            and in what order. Defaults to ``['A', 'B', 'C']``.
        sort_mode (str, optional): Secondary sort criterion. One of
            ``SORT_MODES``. Defaults to ``'valid_days'``.

    Returns:
        dict: A result dictionary with the following keys:

            - ``science_target`` (str)   ``band`` (int)  ``contrast`` (str) ``grade_column`` (str)
            - ``allowed_grades`` (list of str)  ``sort_mode`` (str) ``effective_sort`` (str)  ``sci_mag`` (float or None)
            - ``visibility_pct`` (float)  ``sort_method`` (str)
            - ``observable_windows`` (list of dict), each containing:
                - ``start`` (str)
                - ``end`` (str)
                - ``duration_days`` (float)
                - ``valid_refs`` (list of dict), sorted grade-then-metric, each
                  containing:
                    - ``reference_star`` (str)
                    - ``grade`` (str)
                    - ``grade_rank`` (int)
                    - ``mag`` (float or None)
                    - ``mag_diff`` (float or None): only populated when
                      ``effective_sort='closest_mag'`` and ``sci_mag`` is known
                    - ``n_valid_days`` (int)
                    - ``min_pitch_diff`` (float)
                    - ``valid_dates`` (list of str)

                - ``best_ref`` (dict or None): ``valid_refs[0]`` or None
                - ``pitch_df`` (pandas.DataFrame)
                - ``pitch_csv`` (str): path to the saved CSV file
                - ``avail_df`` (pandas.DataFrame)

            On error, the dict contains an ``'error'`` key with a message
            string instead of ``'observable_windows'``.

    Raises:
        ValueError: If ``(band, contrast)`` is not a recognised combination,
            or if ``sort_mode`` is not a recognised value, or if no valid
            grades remain after filtering.
    """
    key = (band, contrast.lower())
    if key not in GRADE_COLUMNS:
        valid = ', '.join(f"band={b} contrast={ct}" for b, ct in GRADE_COLUMNS)
        raise ValueError(
            f"Unknown (band={band}, contrast='{contrast}'). "
            f"Valid combinations: {valid}"
        )

    if sort_mode not in SORT_MODES:
        raise ValueError(
            f"Unknown sort_mode='{sort_mode}'. Valid: {SORT_MODES}"
        )

    grade_col = GRADE_COLUMNS[key]
    mag_col = 'mag_v' if band == 1 else 'mag_i'
    band_label = 'V' if band == 1 else 'I'

    print(f"\nUsing grade column: {grade_col}  |  mag column: {mag_col}")

    candidates = catalog.copy()
    if 'grade' in candidates.columns:
        candidates['grade'] = candidates['grade'].str.strip()
    elif grade_col in candidates.columns:
        candidates['grade'] = candidates[grade_col].str.strip()
    else:
        raise ValueError(
            f"Neither 'grade' nor '{grade_col}' found in catalog. "
            f"Available: {list(catalog.columns)}"
        )

    active_grades = list(allowed_grades) if allowed_grades else list(REF_GRADES)
    active_grades = [g for g in active_grades if g in REF_GRADES]
    if not active_grades:
        raise ValueError(
            f"allowed_grades={allowed_grades!r} has no valid values. "
            f"Choose from {REF_GRADES}."
        )

    grade_rank_map = {g: i for i, g in enumerate(active_grades)}

    candidates = candidates[candidates['grade'].isin(active_grades)].copy()
    candidates = candidates.dropna(subset=[mag_col])
    candidates['grade_rank'] = (
        candidates['grade']
        .map(grade_rank_map)
        .fillna(99)
        .astype(int)
    )

    print(f"Grade filter: {active_grades} — {len(candidates)} candidate(s) remaining.")

    print(f"Querying SIMBAD for science target '{sci_name}'...")
    coords = get_target_coords([sci_name])
    if sci_name not in coords:
        return {'error': f"Science target '{sci_name}' not found in SIMBAD."}
    sci_coord = coords[sci_name]
    print(f"  Found '{sci_name}'.\n")

    print(f"Looking up {band_label}-band magnitude for '{sci_name}'...")
    sci_mag = get_science_mag(sci_name, band, engine=engine)

    if sort_mode == 'closest_mag' and sci_mag is None:
        effective_sort = 'brightest'
        sort_method = (
            f"grade ({'>'.join(active_grades)}) then brightest {band_label} "
            f"(closest_mag requested but no science target magnitude found)"
        )
    else:
        effective_sort = sort_mode
        sort_method = (
            f"grade ({'>'.join(active_grades)}) then "
            f"{SORT_MODE_LABELS[sort_mode].lower()}"
        )

    print(
        f"\nComputing Roman visibility for '{sci_name}' "
        f"over {analysis_days:.0f} days..."
    )
    times, keepout, _ = compute_keepout(
        {sci_name: sci_coord}, analysis_start, analysis_days, time_step
    )
    sci_keepout = keepout[sci_name]
    visibility_pct = (np.sum(sci_keepout) / len(sci_keepout)) * 100
    print(f"  Observable {visibility_pct:.1f}% of the time.")

    windows = get_observable_windows(times, sci_keepout)
    if not windows:
        return {
            'science_target': sci_name,
            'band': band,
            'contrast': contrast,
            'grade_column': grade_col,
            'allowed_grades': active_grades,
            'sort_mode': sort_mode,
            'effective_sort': effective_sort,
            'sci_mag': sci_mag,
            'error': (
                f"'{sci_name}' is never observable by Roman "
                f"during this period."
            ),
        }

    print(f"\nFound {len(windows)} observable window(s):")
    for i, (_, _, ws, we, wd) in enumerate(windows):
        print(f"  Window {i + 1}: {ws} to {we} ({wd:.1f} days)")

    _, _, _, sci_pitch_full = compute_roman_angles(
        sci_coord, analysis_start, analysis_days, time_step
    )
    sci_pitch_vals = sci_pitch_full.to(u.degree).value

    print("Building reference star coordinates (J2000)...")
    ref_coords = {}
    for _, ref in candidates.iterrows():
        name = ref['main_id']
        if not isinstance(name, str) or name.strip() in SKIP_NAMES:
            continue
        try:
            ref_coords[name] = build_skycoord(ref)
        except Exception as exc:
            print(f"  Warning: could not build coord for '{name}': {exc}")
    print(f"  Successfully built coordinates for {len(ref_coords)} stars.")

    sort_key = make_sort_key(effective_sort)

    results = []
    print(f"\nSort: {sort_method}")
    print("\nSearching for reference stars in each window...")

    for win_idx, (win_start, win_end, ws, we, wd) in enumerate(windows):
        print(f"\nWindow {win_idx + 1}: {ws} to {we} ({wd:.1f} days)")

        win_start_idx = int((win_start.mjd - times[0].mjd) / time_step)
        win_end_idx = int((win_end.mjd - times[0].mjd) / time_step)
        sci_pitch_win = sci_pitch_vals[win_start_idx:win_end_idx + 1]

        n_days_win = int(round(win_end.mjd - win_start.mjd)) + 1
        dates = [
            (win_start + i * u.day).to_value('iso', subfmt='date')
            for i in range(n_days_win)
        ]

        valid_refs = []
        pitch_series = {}

        for _, ref in candidates.iterrows():
            ref_name = ref['main_id']
            if ref_name not in ref_coords:
                continue

            passes, n_days, min_pitch, pd_series, valid_mask = check_ref_in_window(
                ref_coords[ref_name], win_start, win_end, sci_pitch_win,
            )
            pitch_series[ref_name] = pd_series

            if not passes:
                continue

            valid_date_strs = [
                (win_start + int(d) * u.day).to_value('iso', subfmt='date')
                for d in np.where(valid_mask)[0]
            ]

            ref_mag = safe_float(ref.get(mag_col))

            if (
                effective_sort == 'closest_mag'
                and sci_mag is not None
                and ref_mag is not None
            ):
                mag_diff = abs(ref_mag - sci_mag)
            else:
                mag_diff = None

            valid_refs.append({
                'reference_star': ref_name,
                'grade':          ref['grade'],
                'grade_rank':     int(ref['grade_rank']),
                'mag':            ref_mag,
                'mag_diff':       mag_diff,
                'n_valid_days':   n_days,
                'min_pitch_diff': min_pitch,
                'valid_dates':    valid_date_strs,
            })

        valid_refs.sort(key=sort_key)

        avail_data = {
            r['reference_star']: [d in set(r['valid_dates']) for d in dates]
            for r in valid_refs
        }
        avail_df = pd.DataFrame(avail_data, index=dates)
        avail_df.index.name = 'date'

        pitch_data = {}
        for r in valid_refs:
            name = r['reference_star']
            series = pitch_series.get(name, np.array([]))
            vals = list(series[:len(dates)])
            vals += [np.nan] * (len(dates) - len(vals))
            pitch_data[name] = vals

        pitch_df = pd.DataFrame(pitch_data, index=dates)
        pitch_df.index.name = 'date'

        safe_name = sci_name.replace(' ', '_').replace('*', '').strip('_')
        csv_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            f"pitch_table_{safe_name}_band{band}_{contrast}_window{win_idx + 1}.csv",
        )
        pitch_df.to_csv(csv_path, float_format='%.3f')
        print(f"  Pitch table saved → {csv_path}")
        print(f"  Found {len(valid_refs)} valid reference star(s).")

        results.append({
            'start':         ws,
            'end':           we,
            'duration_days': wd,
            'valid_refs':    valid_refs,
            'best_ref':      valid_refs[0] if valid_refs else None,
            'pitch_df':      pitch_df,
            'pitch_csv':     csv_path,
            'avail_df':      avail_df,
        })

    return {
        'science_target':     sci_name,
        'band':               band,
        'contrast':           contrast,
        'grade_column':       grade_col,
        'allowed_grades':     active_grades,
        'sort_mode':          sort_mode,
        'effective_sort':     effective_sort,
        'sci_mag':            sci_mag,
        'visibility_pct':     visibility_pct,
        'sort_method':        sort_method,
        'observable_windows': results,
    }


if __name__ == "__main__":
    import corgidb.ingest

    print("ReferenceStarPicker\n")

    SCIENCE_TARGET = "47 Uma"
    BAND = 1
    CONTRAST = 'high'
    ANALYSIS_START = "2027-01-01T00:00:00"
    ANALYSIS_DAYS = 365
    ALLOWED_GRADES = ['A', 'B', 'C']
    SORT_MODE = 'closest_mag'

    eng = corgidb.ingest.gen_engine('plandb_user', 'plandb_scratch')
    catalog = load_catalog(eng)
    print(f"\nCatalog ready: {len(catalog)} reference stars.\n")

    result = select_ref_star(
        SCIENCE_TARGET, ANALYSIS_START, ANALYSIS_DAYS,
        band=BAND, contrast=CONTRAST,
        catalog=catalog, engine=eng,
        allowed_grades=ALLOWED_GRADES,
        sort_mode=SORT_MODE,
    )

    print("\n" + "=" * 60)
    print(
        f"RESULTS: {result.get('science_target')} "
        f"Band {result.get('band')} {result.get('contrast')}"
    )
    print(f"Allowed grades : {result.get('allowed_grades')}")
    print(f"Sort mode      : {result.get('sort_mode')}")
    sci_mag = result.get('sci_mag')
    print(f"Science mag    : {f'{sci_mag:.2f}' if sci_mag else 'N/A'}")

    if 'error' in result:
        print(f"ERROR: {result['error']}")
    else:
        print(f"Observable {result['visibility_pct']:.1f}% | {result['sort_method']}")
        print("=" * 60)
        for i, win in enumerate(result['observable_windows']):
            print(
                f"\nWindow {i + 1}: {win['start']} → {win['end']} "
                f"({win['duration_days']:.1f} days)"
            )
            for ref in win['valid_refs']:
                if ref['mag_diff'] is not None:
                    mag_str = f"Δmag={ref['mag_diff']:.2f}"
                elif ref['mag'] is not None:
                    mag_str = f"mag={ref['mag']:.2f}"
                else:
                    mag_str = "mag=N/A"
                print(
                    f"  {ref['reference_star']:20s} grade={ref['grade']} "
                    f"{mag_str:14s} "
                    f"{ref['n_valid_days']:3d}d  "
                    f"pitch={ref['min_pitch_diff']:.4f}deg"
                )