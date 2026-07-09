""" Fill missing star data and write a CSV for comparison with the DB.
  1. Load all rows from the Stars table that carry a grade (A/B/C) in any
     of the six grade columns.
  2. For every row, fill in whatever values are missing:
       sy_dist  — derived from sy_plx (parallax in mas to parsec) where absent
       sy_vmag  — queried from SIMBAD where absent
       sy_imag  — queried from SIMBAD; if still absent, derived via EXOSIMS
                  synthetic photometry (sy_vmag + spectype = I-band)
  3. Write one CSV containing only the rows where something was filled.
     The CSV uses the exact same column names as the database so it can
     be diffed directly against the Stars table.
Nothing is written back to the database.
"""

import copy
import json
import os
import importlib.resources
import numpy as np
from datetime import datetime
import pandas as pd
import sqlalchemy as sql
import EXOSIMS.Prototypes.TargetList
from astroquery.simbad import Simbad
from astropy.coordinates import Distance
import astropy.units as u
from synphot import Observation
from synphot.units import VEGAMAG

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

# Persistent on-disk cache — avoids re-querying SIMBAD / re-running synphot
MAG_CACHE_FILE = os.path.join(os.path.dirname(__file__), "mag_cache.json")

# Cache of live EXOSIMS TargetList instances keyed by (mag_v, spec)
_TL_CACHE: dict = {}

# Columns written to the output CSV, in the same order/names as the DB
CSV_COLUMNS = [
    'main_id',
    'st_name',
    'ra',
    'dec',
    'sy_vmag',
    'sy_imag',
    'sy_dist',
    'sy_plx',
    'sy_pmra',
    'sy_pmdec',
    'st_radv',
    'spectype',
    *ALL_GRADE_COLUMNS,
]
# On-disk magnitude cache helpers
def _load_mag_cache() -> dict:
    if not os.path.exists(MAG_CACHE_FILE):
        return {}
    try:
        with open(MAG_CACHE_FILE, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"  Warning: could not read mag cache ({e}). Starting fresh.")
        return {}


def _save_mag_cache(cache: dict) -> None:
    try:
        with open(MAG_CACHE_FILE, 'w') as f:
            json.dump(cache, f, indent=2)
    except Exception as e:
        print(f"  Warning: could not save mag cache ({e}).")

# EXOSIMS synthetic photometry helpers
def _get_targetlist(mag_v: float, spec: str) -> EXOSIMS.Prototypes.TargetList.TargetList:
    """Return (cached) EXOSIMS TargetList for one prototype star."""
    key = (mag_v, spec)
    if key not in _TL_CACHE:
        scriptfile = os.path.join(
            importlib.resources.files("corgietc"), "data", "scripts", "CGI_Noise.json"
        )
        with open(scriptfile, "r") as f:
            specs = json.loads(f.read())
        specs = copy.deepcopy(specs)
        specs["VmagFill"]     = mag_v
        specs["StarSpecFill"] = spec
        _TL_CACHE[key] = EXOSIMS.Prototypes.TargetList.TargetList(**specs)
    return _TL_CACHE[key]


def _synphot_imag(mag_v: float, spec: str) -> float:
    """Derive I-band magnitude from V magnitude + spectral type."""
    TL = _get_targetlist(mag_v, spec)
    template = TL.get_spectral_template(sInd=0, mode=None, Vband=True)
    obs = Observation(template, TL.standard_bands['I'], force='taper')
    return float(obs.effstim(VEGAMAG, vegaspec=TL.OpticalSystem.vega_spectrum).value)

def _simbad_query(names: list) -> dict:
    """Query SIMBAD for sy_vmag (V) and sy_imag (I) for a list of star names.

    Returns:
        dict: {name: {'sy_vmag': float|nan, 'sy_imag': float|nan}}
    """
    if not names:
        return {}

    simbad = Simbad()
    simbad.add_votable_fields('V', 'I')

    def _safe(row, col):
        try:
            v = float(row[col])
            return np.nan if np.isnan(v) else v
        except Exception:
            return np.nan

    out = {}
    for name in names:
        try:
            tbl = simbad.query_object(name)
            if tbl is None or len(tbl) == 0:
                print(f"    SIMBAD: '{name}' not found.")
                continue
            row  = tbl[0]
            vmag = _safe(row, 'V')
            imag = _safe(row, 'I')
            out[name] = {'sy_vmag': vmag, 'sy_imag': imag}
            if not np.isnan(vmag):
                print(f"    {name}: sy_vmag={vmag:.2f} (SIMBAD)")
            if not np.isnan(imag):
                print(f"    {name}: sy_imag={imag:.2f} (SIMBAD)")
        except Exception as e:
            print(f"    SIMBAD query failed for '{name}': {e}")

    return out


def _valid_names(row) -> list[str]:
    """Return ALL usable identifiers for a star row (main_id and st_name).

    Both are treated as independent identifiers — neither is a fallback for
    the other. SIMBAD and cache lookups will be attempted with each one.
    """
    names = []
    for col in ('main_id', 'st_name'):
        val = row.get(col, None)
        if isinstance(val, str) and val.strip() and val.strip() not in SKIP_NAMES:
            names.append(val.strip())
    return names


def _fill_missing(df: pd.DataFrame, was_missing: dict) -> pd.DataFrame:
    """Fill sy_dist, sy_vmag, sy_imag in df where they are NaN.

    was_missing is a dict of {col: boolean Series} capturing which cells
    were NULL in the DB *before* any filling — passed in from load_and_fill
    so the CSV diff is accurate even when values come from the cache.

    Fill order per column:
      sy_dist  — parallax conversion (done before this call)
      sy_vmag  — cache → SIMBAD (tries main_id AND st_name independently)
      sy_imag  — cache → SIMBAD (both names) → synthetic photometry

    Modifies df in-place and returns it.
    """
    cache         = _load_mag_cache()
    cache_updated = False

    # sy_vmag / sy_imag: on-disk cache 
    # Try BOTH main_id and st_name as independent cache keys — merge results
    n_cache = 0
    for idx, row in df[df['sy_vmag'].isna() | df['sy_imag'].isna()].iterrows():
        for name in _valid_names(row):
            entry = cache.get(name, {})
            if pd.isna(df.at[idx, 'sy_vmag']) and 'sy_vmag' in entry:
                df.at[idx, 'sy_vmag'] = entry['sy_vmag']
                n_cache += 1
            if pd.isna(df.at[idx, 'sy_imag']) and 'sy_imag' in entry:
                df.at[idx, 'sy_imag'] = entry['sy_imag']
                n_cache += 1
    if n_cache:
        print(f"  Magnitudes: restored {n_cache} value(s) from cache.")

    # sy_vmag / sy_imag: SIMBAD 
    # Collect ALL valid identifiers (main_id + st_name) for every row still
    # missing data. Each identifier is queried independently — results from
    # either name are merged back into the same row.
    name_to_indices: dict[str, list[int]] = {}
    for idx, row in df[df['sy_vmag'].isna() | df['sy_imag'].isna()].iterrows():
        for name in _valid_names(row):
            name_to_indices.setdefault(name, []).append(idx)

    if name_to_indices:
        unique_names = list(name_to_indices.keys())
        print(f"  Querying SIMBAD for {len(unique_names)} identifier(s) "
              f"(main_id + st_name combined, {len(df[df['sy_vmag'].isna() | df['sy_imag'].isna()])} rows)...")
        simbad_results = _simbad_query(unique_names)
        for name, mags in simbad_results.items():
            for idx in name_to_indices.get(name, []):
                entry = cache.setdefault(name, {})
                if pd.isna(df.at[idx, 'sy_vmag']) and not np.isnan(mags['sy_vmag']):
                    df.at[idx, 'sy_vmag'] = mags['sy_vmag']
                    entry['sy_vmag'] = mags['sy_vmag']
                    cache_updated    = True
                if pd.isna(df.at[idx, 'sy_imag']) and not np.isnan(mags['sy_imag']):
                    df.at[idx, 'sy_imag'] = mags['sy_imag']
                    entry['sy_imag'] = mags['sy_imag']
                    cache_updated    = True

    # sy_imag: synthetic photometry (all remaining, including those whose sy_vmag was just filled by SIMBAD above) 
    still = int(df['sy_imag'].isna().sum())
    if still:
        candidates = df[df['sy_imag'].isna() & df['sy_vmag'].notna()]
        print(f"  sy_imag: {still} still missing, "
              f"{len(candidates)} have sy_vmag available for synphot...")

        # Diagnose how many will be skipped due to missing spectype
        n_no_spec = int((
            candidates['spectype'].isna() |
            candidates['spectype'].astype(str).str.strip().isin(['', 'nan'])
        ).sum())
        print(f"    {len(candidates)} candidates: "
              f"{len(candidates) - n_no_spec} have spectype, "
              f"{n_no_spec} missing spectype (will be skipped).")
        if n_no_spec > 0:
            missing_spec = candidates[
                candidates['spectype'].isna() |
                candidates['spectype'].astype(str).str.strip().isin(['', 'nan'])
            ]
            print("    Sample of stars skipped (no spectype):")
            print(missing_spec[['main_id', 'st_name', 'sy_vmag', 'spectype']]
                  .head(10).to_string(index=False))

        n_ok = n_fail = 0
        for idx, row in candidates.iterrows():
            spec = row.get('spectype', None)
            if not isinstance(spec, str) or not spec.strip() or spec.strip() == 'nan':
                n_fail += 1
                continue
            try:
                imag = _synphot_imag(float(row['sy_vmag']), spec.strip())
                df.at[idx, 'sy_imag'] = imag
                # Cache under whichever name is available
                # Cache the result under every valid identifier for this row
                for n in _valid_names(row):
                    cache.setdefault(n, {})['sy_imag'] = imag
                    cache_updated = True
                label = ' / '.join(_valid_names(row)) or '(unknown)'
                n_ok += 1
                print(f"    {label}: sy_imag={imag:.2f} "
                      f"(sy_vmag={row['sy_vmag']:.2f}, spec={spec.strip()})")
            except Exception as e:
                label = ' / '.join(_valid_names(row)) or '(unknown)'
                print(f"    {label}: synphot failed ({e})")
                n_fail += 1

        still_no_vmag = int((df['sy_imag'].isna() & df['sy_vmag'].isna()).sum())
        print(f"  Synthetic photometry: {n_ok} succeeded, {n_fail} skipped/failed.")
        if still_no_vmag:
            print(f"  {still_no_vmag} star(s) still missing sy_imag because "
                  f"sy_vmag is also missing (no V-band reference available).")

    if cache_updated:
        _save_mag_cache(cache)
        print(f"  Cache saved → {MAG_CACHE_FILE}")

    return df


def load_and_fill(
    engine: sql.engine.base.Engine,
    output_csv: str | None = None,
) -> pd.DataFrame:
    """Load all graded stars, fill missing values, write gap-fill CSV.

    The output CSV contains only the rows where at least one value was filled,
    and only shows the filled cells — cells already present in the DB are left
    blank so the file reads as a pure diff against the database.

    Args:
        engine:     SQLAlchemy engine connected to plandb.
        output_csv: Destination path for the CSV.
                    Defaults to filled_catalog.csv next to this file.

    Returns:
        pd.DataFrame — the complete filled catalog (all rows, all columns),
        using the same column names as the DB.  Suitable for use by picker.py.
    """
    if output_csv is None:
        timestamp  = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_csv = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            f"filled_catalog_{timestamp}.csv"
        )

    # Load from DB
    metadata    = sql.MetaData()
    stars_table = sql.Table('Stars', metadata, autoload_with=engine)

    # Use the actual column names reflected from the DB (hasattr is unreliable
    # with SQLAlchemy — it returns True for any attribute, not just real columns)
    actual_cols = {c.name for c in stars_table.columns}
    print(f"  Stars table columns: {sorted(actual_cols)}")

    present_grade_cols = [col for col in ALL_GRADE_COLUMNS if col in actual_cols]
    absent_grade_cols  = set(ALL_GRADE_COLUMNS) - set(present_grade_cols)
    if absent_grade_cols:
        print(f"  Note: grade column(s) not yet in DB: {absent_grade_cols}")

    # If ALL six new grade columns are fully NULL, fall back to st_psfgrade
    all_null = True
    if present_grade_cols:
        with engine.connect() as _conn:
            for _col in present_grade_cols:
                _row = _conn.execute(
                    sql.text(f"SELECT COUNT(*) FROM Stars WHERE {_col} IS NOT NULL")
                ).scalar()
                if _row and _row > 0:
                    all_null = False
                    break

    if not present_grade_cols or all_null:
        fallback = 'st_psfgrade'
        if fallback in actual_cols:
            print(f"  New grade columns are all NULL — falling back to '{fallback}'.")
            present_grade_cols = [fallback]
        else:
            raise RuntimeError(
                f"No populated grade columns found. "
                f"Tried {ALL_GRADE_COLUMNS} (all NULL) and '{fallback}' (not found). "
                f"Available columns: {sorted(actual_cols)}"
            )

    # Load ALL stars — we want to fill missing sy_imag and sy_dist
    # for every row in the table, not just graded ones.
    # We still select the grade columns so they appear in the CSV output.
    grade_sel_cols = [
        getattr(stars_table.c, col)
        for col in present_grade_cols
        if col in actual_cols
    ]
    stmt = sql.select(
        stars_table.c.main_id,
        stars_table.c.st_name,
        stars_table.c.ra,
        stars_table.c.dec,
        stars_table.c.sy_vmag,
        stars_table.c.sy_imag,
        stars_table.c.sy_dist,
        stars_table.c.sy_plx,
        stars_table.c.sy_pmra,
        stars_table.c.sy_pmdec,
        stars_table.c.st_radv,
        stars_table.c.spectype,
        *grade_sel_cols,
    )  # NO WHERE clause — load every star

    with engine.connect() as conn:
        df = pd.read_sql(stmt, conn)

    # Add placeholder columns for any grade columns not in this DB
    for col in ALL_GRADE_COLUMNS:
        if col not in df.columns:
            df[col] = np.nan

    # Coerce numeric types
    for col in ('sy_vmag', 'sy_imag', 'sy_dist', 'sy_plx',
                'sy_pmra', 'sy_pmdec', 'st_radv'):
        df[col] = pd.to_numeric(df[col], errors='coerce')

    print(f"Loaded {len(df)} star(s) from Stars table (all rows).")
    print(f"  Missing sy_vmag : {df['sy_vmag'].isna().sum()}")
    print(f"  Missing sy_imag : {df['sy_imag'].isna().sum()}")
    print(f"  Missing sy_dist : {df['sy_dist'].isna().sum()}")
    print(f"  Has sy_plx      : {df['sy_plx'].notna().sum()}")
    # If sy_dist is missing but sy_plx is available, derive it now
    # before the snapshot so it appears correctly in the filled count
    # Note: we derive dist HERE so was_missing captures the true DB state
    pre_dist_missing = df['sy_dist'].isna().copy()
    dist_derivable = pre_dist_missing & df['sy_plx'].notna() & (df['sy_plx'] > 0)
    n_dist = int(dist_derivable.sum())
    if n_dist:
        df.loc[dist_derivable, 'sy_dist'] = Distance(
            parallax=df.loc[dist_derivable, 'sy_plx'].values * u.mas
        ).pc
        print(f"  Derived sy_dist from sy_plx for {n_dist} star(s).")

    # Which cells are NULL in the DB, BEFORE any filling 
    # sy_dist is already derived above; use pre_dist_missing to record
    # what was truly absent in the DB so the CSV diff is accurate.
    was_missing = {
        'sy_vmag': df['sy_vmag'].isna().copy(),
        'sy_imag': df['sy_imag'].isna().copy(),
        'sy_dist': pre_dist_missing,
    }
    df = _fill_missing(df, was_missing)

    print(f"\nAfter filling:")
    print(f"  Still missing sy_vmag : {df['sy_vmag'].isna().sum()}")
    print(f"  Still missing sy_imag : {df['sy_imag'].isna().sum()}")
    print(f"  Still missing sy_dist : {df['sy_dist'].isna().sum()}")

    # Include every row where at least one value was missing in the DB and
    # has now been filled.  Show ALL columns so each row has full context
    # (st_name, coordinates, grade columns, etc.).  The filled values
    # (sy_vmag, sy_imag, sy_dist) appear in their normal columns.
    # An extra 'filled_columns' column lists which fields were gap-filled
    # for that row, making it easy to see at a glance what was added.
    any_filled = pd.Series(False, index=df.index)
    for col in ('sy_vmag', 'sy_imag', 'sy_dist'):
        any_filled |= was_missing[col] & df[col].notna()

    out = df[any_filled].copy()

    # Add a 'filled_columns' column listing every field that was NULL in the
    # DB and has now been given a value — e.g. 'sy_vmag, sy_imag, sy_dist'
    fill_cols = ('sy_vmag', 'sy_imag', 'sy_dist')
    labels = []
    for idx in out.index:
        filled = [
            col for col in fill_cols
            if was_missing[col].loc[idx] and pd.notna(out.at[idx, col])
        ]
        labels.append(', '.join(filled) if filled else 'none')
    out['filled_columns'] = labels

    # Also add individual boolean flag columns so it's easy to filter in Excel
    for col in fill_cols:
        out[f'{col}_filled'] = [
            was_missing[col].loc[idx] and pd.notna(out.at[idx, col])
            for idx in out.index
        ]

    # Column order: identity columns first, then filled_columns indicator,
    # then the three fillable values, then everything else
    priority = ['main_id', 'st_name', 'filled_columns',
                'sy_vmag_filled', 'sy_imag_filled', 'sy_dist_filled',
                'sy_vmag', 'sy_imag', 'sy_dist', 'sy_plx',
                'ra', 'dec', 'sy_pmra', 'sy_pmdec', 'st_radv', 'spectype',
                *ALL_GRADE_COLUMNS]
    col_order = [c for c in priority if c in out.columns]
    col_order += [c for c in out.columns if c not in col_order]
    out = out[col_order]

    out.to_csv(output_csv, index=False, float_format='%.6f')
    print(f"\nGap-fill CSV: {len(out)} row(s) written → {output_csv}")
    print(f"  Columns filled: sy_vmag={was_missing['sy_vmag'].sum()} missing, "
          f"sy_imag={was_missing['sy_imag'].sum()} missing, "
          f"sy_dist={was_missing['sy_dist'].sum()} missing in DB.")

    return df


if __name__ == "__main__":
    import corgidb.ingest

    eng = corgidb.ingest.gen_engine('plandb_user', 'plandb_scratch')

    load_and_fill(eng)
