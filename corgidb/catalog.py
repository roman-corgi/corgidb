"""
Reference star catalog loading and magnitude resolution for ReferenceStarPicker.
  - Loading from plandb_scratch Stars table
  - Querying SIMBAD for missing V and I magnitudes
  - Deriving missing I-band magnitudes from spectral type via synthetic photometry
  - Querying SIMBAD for the science target's magnitude

MAGNITUDE LOOKUP CHAIN
  1. plandb Stars table  (sy_vmag / sy_imag)
  2. SIMBAD query        (V / I flux columns)
  3. Synthetic photometry: instantiates a real EXOSIMS TargetList via the
     CGI_Noise.json spec with VmagFill set to the known V-band magnitude and
     StarSpecFill set to the spectral type, then calls get_spectral_template
     and integrates through the Cousins-I filter.

DISTANCE
  Missing sy_dist values are derived from sy_plx (parallax in mas) using
  astropy.coordinates.Distance for the exact conversion.
"""

import copy
import json
import os
import importlib.resources
import numpy as np
import pandas as pd
import sqlalchemy as sql
import EXOSIMS.Prototypes.TargetList
from astroquery.simbad import Simbad
from astropy.coordinates import Distance
import astropy.units as u
from synphot import Observation
from synphot.units import VEGAMAG

# Grades to treat as reference star candidates
REF_GRADES = ['A', 'B', 'C']
SKIP_NAMES = {'-', 'TBD', '?', ''}

# Persistent on-disk cache for resolved magnitudes (SIMBAD + synphot results).
# Only entries that were *not* present in the database are stored here, so the
# file stays small and the database always wins on the next load.
MAG_CACHE_FILE = os.path.join(os.path.dirname(__file__), "mag_cache.json")

# Cache of instantiated TargetLists keyed by (mag_v, spec)
_TL_CACHE: dict = {}



def _load_mag_cache() -> dict:
    """Load the on-disk magnitude cache.

    Returns:
        dict keyed by star name; values are dicts with optional keys
        'mag_v' and 'mag_i' (both float).  Returns {} if the file does
        not exist or cannot be parsed.
    """
    if not os.path.exists(MAG_CACHE_FILE):
        return {}
    try:
        with open(MAG_CACHE_FILE, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"  Warning: could not read mag cache ({e}). Starting fresh.")
        return {}


def _save_mag_cache(cache: dict) -> None:
    """Write the magnitude cache to disk.

    Args:
        cache (dict): The full cache dict to serialise.
    """
    try:
        with open(MAG_CACHE_FILE, 'w') as f:
            json.dump(cache, f, indent=2)
    except Exception as e:
        print(f"  Warning: could not save mag cache ({e}).")


def _get_targetlist(mag_v: float, spec: str) -> EXOSIMS.Prototypes.TargetList.TargetList:
    """Return a real TargetList for a single prototype star.

    Uses the CGI_Noise.json default spec with VmagFill and StarSpecFill
    overridden to represent the star of interest. Results are cached so
    that repeated calls with the same (mag_v, spec) do not re-instantiate.

    Args:
        mag_v (float): V-band magnitude of the star.
        spec (str):    Normalised spectral type string e.g. 'G2V'.

    Returns:
        EXOSIMS TargetList with one star (sInd=0).
    """
    key = (mag_v, spec)
    if key not in _TL_CACHE:
        scriptfile = os.path.join(
            importlib.resources.files("corgietc"), "data", "scripts", "CGI_Noise.json"
        )
        with open(scriptfile, "r") as f:
            specs = json.loads(f.read())

        specs = copy.deepcopy(specs)
        specs["VmagFill"] = mag_v
        specs["StarSpecFill"] = spec

        _TL_CACHE[key] = EXOSIMS.Prototypes.TargetList.TargetList(**specs)

    return _TL_CACHE[key]


def _get_imag_synphot(mag_v: float, spec: str) -> float:
    """Derive I-band magnitude for one star using EXOSIMS synthetic photometry.

    Args:
        mag_v (float): Known V-band magnitude.
        spec (str):    Spectral type string e.g. 'G2V'.

    Returns:
        float: Derived I-band magnitude.
    """
    TL = _get_targetlist(mag_v, spec)
    template_renorm = TL.get_spectral_template(sInd=0, mode=None, Vband=True)
    obs = Observation(template_renorm, TL.standard_bands['I'], force='taper')
    return float(
        obs.effstim(VEGAMAG, vegaspec=TL.OpticalSystem.vega_spectrum).value
    )


def _simbad_query_mags(names: list) -> dict:
    """Query SIMBAD for V and I(Cousins) magnitudes for a list of star names.

    Args:
        names (list[str]): Star names to query.

    Returns:
        dict keyed by star name; values are dicts with 'mag_v' and 'mag_i'
        (float or NaN). Stars not found in SIMBAD are omitted.
    """
    if not names:
        return {}

    simbad = Simbad()
    simbad.add_votable_fields('V', 'I')

    def safe(row, col):
        try:
            v = float(row[col])
            return np.nan if np.isnan(v) else v
        except Exception:
            return np.nan

    results = {}
    for name in names:
        try:
            tbl = simbad.query_object(name)
            if tbl is None or len(tbl) == 0:
                print(f"    SIMBAD: '{name}' not found.")
                continue

            row = tbl[0]
            mag_v = safe(row, 'V')
            mag_i = safe(row, 'I')
            results[name] = {'mag_v': mag_v, 'mag_i': mag_i}

            if not np.isnan(mag_v):
                print(f"    {name}: V={mag_v:.2f} from SIMBAD")
            if not np.isnan(mag_i):
                print(f"    {name}: I={mag_i:.2f} from SIMBAD")

        except Exception as e:
            print(f"    SIMBAD query failed for '{name}': {e}")

    return results


def _derive_imag_from_synphot(df: pd.DataFrame) -> pd.DataFrame:
    """Derive missing I-band magnitudes using EXOSIMS synthetic photometry.

    Args:
        df (pd.DataFrame):
            Catalog with columns: 'name', 'mag_v', 'mag_i', 'spectype'.

    Returns:
        pd.DataFrame with mag_i filled where possible.
    """
    needs_i = df[df['mag_i'].isna() & df['mag_v'].notna()].copy()
    if needs_i.empty:
        return df

    print(f"  Deriving I-band via synthetic photometry for {len(needs_i)} stars...")

    n_derived = 0
    n_failed  = 0

    for idx, row in needs_i.iterrows():
        spec = row.get('spectype', None)
        if not isinstance(spec, str) or not spec.strip():
            n_failed += 1
            continue

        try:
            spec_norm = spec.strip()
            mag_i = _get_imag_synphot(row['mag_v'], spec_norm)
            df.at[idx, 'mag_i'] = mag_i
            n_derived += 1
            print(f"    {row['name']}: I={mag_i:.2f} "
                  f"(V={row['mag_v']:.2f}, spec={spec_norm})")

        except Exception as e:
            print(f"    {row['name']}: synthetic photometry failed ({e})")
            n_failed += 1

    print(f"  Synthetic photometry: {n_derived} succeeded, {n_failed} failed.")
    return df


def _resolve_magnitudes(df: pd.DataFrame) -> pd.DataFrame:
    """Fill missing V and I magnitudes using a three-step fallback chain.

    Step 0 - Persistent on-disk cache (mag_cache.json) for values resolved
             in previous runs via SIMBAD or synthetic photometry.
    Step 1 - Database values already in df (mag_v / mag_i from plandb).
             Database values always overwrite the cache for the same star.
    Step 2 - SIMBAD query for any star still missing mag_v or mag_i.
    Step 3 - Synthetic photometry via EXOSIMS TargetList.get_spectral_template.

    New values found in Steps 2 and 3 are written back to the cache so
    subsequent runs skip the expensive external lookups for known stars.

    Args:
        df (pd.DataFrame):
            Catalog with columns: 'name', 'mag_v', 'mag_i', 'spectype'.

    Returns:
        pd.DataFrame with mag_v and mag_i filled as completely as possible.
    """
    mag_cache = _load_mag_cache()
    cache_updated = False

    # Step 0: apply cache for stars still missing mag_v or mag_i
    needs_fill = df['mag_v'].isna() | df['mag_i'].isna()
    n_from_cache = 0
    for idx, row in df[needs_fill].iterrows():
        name = row['name']
        if not isinstance(name, str) or name.strip() in SKIP_NAMES:
            continue
        entry = mag_cache.get(name, {})
        if df.at[idx, 'mag_v'] != df.at[idx, 'mag_v'] and 'mag_v' in entry:  # isnan check
            df.at[idx, 'mag_v'] = entry['mag_v']
            n_from_cache += 1
        if df.at[idx, 'mag_i'] != df.at[idx, 'mag_i'] and 'mag_i' in entry:
            df.at[idx, 'mag_i'] = entry['mag_i']
            n_from_cache += 1

    if n_from_cache:
        print(f"  Restored {n_from_cache} magnitude value(s) from cache "
              f"({MAG_CACHE_FILE}).")

    # Step 2: SIMBAD — only for stars still missing after cache
    needs_simbad = [
        n for n in df.loc[
            df['mag_v'].isna() | df['mag_i'].isna(), 'name'
        ].tolist()
        if isinstance(n, str) and n.strip() not in SKIP_NAMES
    ]

    if needs_simbad:
        print(f"  Querying SIMBAD for {len(needs_simbad)} stars "
              f"missing V or I mag...")
        simbad_data = _simbad_query_mags(needs_simbad)

        for name, mags in simbad_data.items():
            mask = df['name'] == name
            entry = mag_cache.setdefault(name, {})
            if df.loc[mask, 'mag_v'].isna().any() and not np.isnan(mags['mag_v']):
                df.loc[mask, 'mag_v'] = mags['mag_v']
                entry['mag_v'] = mags['mag_v']
                cache_updated = True
            if df.loc[mask, 'mag_i'].isna().any() and not np.isnan(mags['mag_i']):
                df.loc[mask, 'mag_i'] = mags['mag_i']
                entry['mag_i'] = mags['mag_i']
                cache_updated = True

    # Step 3: Synthetic photometry V -> I — only for stars still missing
    still_missing_i = df['mag_i'].isna().sum()
    if still_missing_i > 0:
        print(f"  {still_missing_i} stars still missing I-band — "
              f"attempting synthetic photometry derivation...")
        # Capture which rows had mag_i before so we can detect new values
        had_i = df['mag_i'].notna().copy()
        df = _derive_imag_from_synphot(df)
        newly_filled = df['mag_i'].notna() & ~had_i
        for idx, row in df[newly_filled].iterrows():
            name = row['name']
            if isinstance(name, str) and name.strip() not in SKIP_NAMES:
                mag_cache.setdefault(name, {})['mag_i'] = float(row['mag_i'])
                cache_updated = True

    if cache_updated:
        _save_mag_cache(mag_cache)
        print(f"  Magnitude cache updated → {MAG_CACHE_FILE}")

    return df


def load_catalog(engine: sql.engine.base.Engine) -> pd.DataFrame:
    """Load the reference star catalog from plandb, filling missing magnitudes.

    Queries Stars where st_psfgrade is A, B, or C. Uses J2000 coordinates
    (ra, dec). Derives sy_dist from sy_plx where sy_dist is NULL using
    astropy.coordinates.Distance.

    Magnitude resolution chain:
      1. Database values (sy_vmag / sy_imag)
      2. SIMBAD query for any missing V or I
      3. Synthetic photometry V->I derivation via EXOSIMS TargetList

    Args:
        engine (sqlalchemy.engine.base.Engine):
            Connected engine for plandb.

    Returns:
        pd.DataFrame with columns:
            name, st_name, grade, mag_v, mag_i, ra, dec,
            sy_dist, sy_pmra, sy_pmdec, st_radv, spectype, sy_plx
    """
    metadata    = sql.MetaData()
    stars_table = sql.Table('Stars', metadata, autoload_with=engine)

    stmt = sql.select(
        stars_table.c.main_id,
        stars_table.c.st_name,
        stars_table.c.st_psfgrade,
        stars_table.c.sy_vmag,
        stars_table.c.sy_imag,
        stars_table.c.ra,
        stars_table.c.dec,
        stars_table.c.sy_pmra,
        stars_table.c.sy_pmdec,
        stars_table.c.sy_dist,
        stars_table.c.sy_plx,
        stars_table.c.st_radv,
        stars_table.c.spectype,
    ).where(stars_table.c.st_psfgrade.in_(REF_GRADES))

    with engine.connect() as conn:
        df = pd.read_sql(stmt, conn)

    df = df.rename(columns={
        'main_id':     'name',
        'st_psfgrade': 'grade',
        'sy_vmag':     'mag_v',
        'sy_imag':     'mag_i',
    })

    for col in ('mag_v', 'mag_i', 'sy_pmra', 'sy_pmdec', 'sy_dist', 'sy_plx',
                'st_radv'):
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Derive distance from parallax where sy_dist is NULL
    missing_dist = (
        df['sy_dist'].isna() & df['sy_plx'].notna() & (df['sy_plx'] > 0)
    )
    n_derived = missing_dist.sum()
    if n_derived:
        df.loc[missing_dist, 'sy_dist'] = Distance(
            parallax=df.loc[missing_dist, 'sy_plx'].values * u.mas
        ).pc
        print(f"  Derived sy_dist from sy_plx for {n_derived} stars.")

    df[['sy_pmra', 'sy_pmdec', 'st_radv']] = (
        df[['sy_pmra', 'sy_pmdec', 'st_radv']].fillna(0.0)
    )
    df['grade'] = df['grade'].str.strip()

    print(f"Loaded {len(df)} reference stars from database (J2000 coordinates).")
    print(f"  Missing V-band : {df['mag_v'].isna().sum()}")
    print(f"  Missing I-band : {df['mag_i'].isna().sum()}")
    print(f"  Missing dist   : {df['sy_dist'].isna().sum()}")

    df = _resolve_magnitudes(df)

    print(f"After resolution — still missing V: {df['mag_v'].isna().sum()} "
          f"| still missing I: {df['mag_i'].isna().sum()}")

    still_missing = df[df['mag_i'].isna()]
    if not still_missing.empty:
        print(f"\nRemaining I-band gaps by spectral type "
              f"({len(still_missing)} stars):")
        print(still_missing['spectype'].value_counts(dropna=False).to_string())

    return df


def get_science_mag(
    sci_name: str,
    band: int,
    engine: sql.engine.base.Engine = None,
) -> float | None:
    """Query the database then SIMBAD for the science target's magnitude.

    Fallback chain:
      1. plandb Stars table (sy_vmag / sy_imag) matched on st_name
         (case-insensitive) — if engine is provided
      2. SIMBAD query (V / I flux columns)
      3. Returns None  caller falls back to brightest-first sorting

    Args:
        sci_name (str):  Science target name recognized by SIMBAD.
        band (int):      1 = V band, 3 or 4 = I band.
        engine:          SQLAlchemy engine for plandb (optional).

    Returns:
        float or None: Magnitude if found; None triggers brightest-first sorting.
    """
    band_label = 'V' if band == 1 else 'I'
    db_col     = 'sy_vmag' if band == 1 else 'sy_imag'
    flux_col   = 'V'       if band == 1 else 'I'

    # Step 1: Database
    if engine is not None:
        try:
            metadata    = sql.MetaData()
            stars_table = sql.Table('Stars', metadata, autoload_with=engine)
            stmt = sql.select(
                getattr(stars_table.c, db_col)
            ).where(
                sql.func.lower(stars_table.c.st_name) == sci_name.lower()
            )

            with engine.connect() as conn:
                row = conn.execute(stmt).fetchone()

            if row is not None:
                val = row[0]
                if val is not None and not (isinstance(val, float) and np.isnan(val)):
                    mag = float(val)
                    print(f"  {band_label}-band magnitude from database: {mag:.2f}")
                    return mag
                print(f"  '{sci_name}' found in database but {db_col} is NULL.")
            else:
                print(f"  '{sci_name}' not found in Stars table.")
        except Exception as e:
            print(f"  Database magnitude lookup failed: {e}")

    # Step 2: SIMBAD
    print(f"  Querying SIMBAD for {band_label}-band magnitude...")
    try:
        simbad = Simbad()
        simbad.add_votable_fields('V', 'I')
        tbl = simbad.query_object(sci_name)

        if tbl is None or len(tbl) == 0:
            print(f"  SIMBAD returned no results for '{sci_name}'.")
        else:
            val = tbl[0][flux_col]
            if hasattr(val, 'mask') and val.mask:
                mag = np.nan
            else:
                mag = float(val)

            if not np.isnan(mag):
                print(f"  {band_label}-band magnitude from SIMBAD: {mag:.2f}")
                return mag
    except Exception as e:
        print(f"  SIMBAD magnitude query failed: {e}")

    # Step 3: Give up
    print(f"  No {band_label} magnitude found. "
          f"Falling back to brightest-first sorting within grade.")
    return None