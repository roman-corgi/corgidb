"""
Reference star selection logic for ReferenceStarPicker.
  - Building astropy SkyCoords from catalog rows using J2000 coordinates
  - Computing Roman's observable windows via keepout
  - Checking solar angle and pitch angle constraints day-by-day per window
  - Returning ALL valid reference stars per window, ranked by observable days

OVERALL FLOW
1. compute_keepout() on the science target
2. For each observable window, check all candidate reference stars ranked by:
       a. Grade (A is best, then B, then C)
       b. Closest magnitude to the science target in the chosen band (V or I).
          If no science target magnitude is available, sort brightest-first.

3. For each candidate, check TWO constraints DAY BY DAY:
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
import corgidb.ingest
from catalog import REF_GRADES, SKIP_NAMES, get_science_mag, load_catalog

from roman_pointing.roman_observability import (
    get_target_coords,
    compute_roman_angles,
    compute_keepout,
)
SUN_MIN        = 54    # minimum solar angle (deg)
SUN_MAX        = 126   # maximum solar angle (deg)
MAX_PITCH_DIFF = 5.0   # maximum pitch angle difference sci <-> ref (deg)


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

    pitch_diff_series              = pitch_diff.copy().astype(float)
    pitch_diff_series[~solar_ok]   = np.nan   # NaN = solar angle blocked

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
    catalog: pd.DataFrame,
    engine=None,
    time_step: float = 1.0,
) -> dict:
    """Find all valid reference stars for each observable window of a science target.

    Reference stars are checked day-by-day within each window. A star is
    valid if it passes both the solar angle and pitch angle constraints on
    at least one day. Results are sorted by number of valid days (descending).

    Args:
        sci_name (str):         Science target name recognized by SIMBAD
        analysis_start (str):   ISO start date, e.g. '2027-01-01T00:00:00'
        analysis_days (float):  Total span to analyze (days)
        band (int):             1 = V band, 3 or 4 = I band
        catalog (pd.DataFrame): Reference star catalog from load_catalog()
        engine:                 SQLAlchemy engine for plandb (optional)
        time_step (float):      Time resolution for angle calculations (days)

    Returns:
        dict with keys:
            'science_target', 'band', 'visibility_pct', 'sort_method',
            'observable_windows' (list of per-window result dicts), each with:
                'start', 'end', 'duration_days',
                'valid_refs' (list of all passing ref stars, sorted by n_valid_days),
                'best_ref'   (ref star with most valid days, or None)
    """
    if band == 1:
        mag_col, band_label = 'mag_v', 'V'
    elif band in (3, 4):
        mag_col, band_label = 'mag_i', 'I'
    else:
        return {'error': f"Unknown band '{band}'. Use 1 (V) or 3/4 (I)."}

    # Science target coords
    print(f"Querying SIMBAD for science target '{sci_name}'...")
    coords = get_target_coords([sci_name])
    if sci_name not in coords:
        return {'error': f"Science target '{sci_name}' not found in SIMBAD."}
    sci_coord = coords[sci_name]
    print(f"  Found '{sci_name}'.\n")

    # Science target magnitude
    print(f"Looking up {band_label}-band magnitude for '{sci_name}'...")
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
    candidates  = catalog.dropna(subset=[mag_col]).copy()
    grade_order = {g: i for i, g in enumerate(REF_GRADES)}
    candidates['grade_rank'] = candidates['grade'].map(grade_order).fillna(99)

    if sci_mag is not None:
        candidates['mag_diff'] = (candidates[mag_col] - sci_mag).abs()
        candidates  = candidates.sort_values(['grade_rank', 'mag_diff'])
        sort_method = (f"grade (A>B>C) then closest {band_label}-band "
                       f"magnitude to science target")
    else:
        candidates['mag_diff'] = np.nan
        candidates  = candidates.sort_values(['grade_rank', mag_col])
        sort_method = f"grade (A>B>C) then brightest {band_label}"

    print(f"\nSorting by: {sort_method}")
    print(f"Candidates available: {len(candidates)}")

    # Build SkyCoords (J2000 only)
    print("Building reference star coordinates (J2000)...")
    ref_coords = {}
    for _, ref in candidates.iterrows():
        name = ref['name']
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
                ref_name = ref['name']
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

        # Save to CSV — never printed, open in Excel or query in pandas
        safe_name = sci_name.replace(' ', '_').replace('*', '').strip('_')
        csv_path  = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 f"pitch_table_{safe_name}_window{win_idx + 1}.csv")
        pitch_df.to_csv(csv_path, float_format='%.3f')
        print(f"  Pitch table saved → {csv_path}")

        print(f"  Found {len(valid_refs)} valid reference star(s) for this window.")
        results.append({
            'start':         ws,
            'end':           we,
            'duration_days': wd,
            'valid_refs':    valid_refs,
            'best_ref':      valid_refs[0] if valid_refs else None,
            'pitch_df':      pitch_df,    # DataFrame available for programmatic use
            'pitch_csv':     csv_path,
        })

    return {
        'science_target':     sci_name,
        'band':               band,
        'visibility_pct':     visibility_pct,
        'sort_method':        sort_method,
        'observable_windows': results,
    }


if __name__ == "__main__":

    print("ReferenceStarPicker\n")

    SCIENCE_TARGET = "14 Her"
    BAND           = 1
    ANALYSIS_START = "2027-01-01T00:00:00"
    ANALYSIS_DAYS  = 365

    eng = corgidb.ingest.gen_engine('plandb_user', 'plandb_scratch')

    catalog = load_catalog(eng)
    print(f"\nCatalog ready: {len(catalog)} reference stars.\n")

    result = select_ref_star(
        SCIENCE_TARGET, ANALYSIS_START, ANALYSIS_DAYS, BAND, catalog,
        engine=eng,
    )

    print("\n" + "=" * 60)
    print(f"RESULTS FOR: {result.get('science_target')} "
          f"(Band {result.get('band')})")

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
                          f"min pitch diff={ref['min_pitch_diff']:.2f} deg)")