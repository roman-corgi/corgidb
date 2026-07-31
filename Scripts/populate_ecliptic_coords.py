"""Populate missing ecliptic coordinates and uncertainties in the Stars table.

For each row missing elat/elon (with ra/dec present), transform the ICRS
coordinates to BarycentricMeanEcliptic. Kinematic columns (sy_plx,
sy_pmra, sy_pmdec, st_radv) are folded in when available.

Uncertainties (elaterr1/elaterr2, elonerr1/elonerr2) are obtained by
Monte Carlo, drawing from a two-piece normal distribution that honors
the asymmetric err1 (upper) and err2 (lower) columns on every input
axis. Rows whose ra/dec uncertainties are all NULL still get elat/elon
populated, but their ecliptic uncertainty columns remain NULL.

The script is dry-run by default; pass --write to commit.

Example:
    python -m Scripts.populate_ecliptic_coords --username plandb_user
    python -m Scripts.populate_ecliptic_coords --username plandb_user --write
"""

import argparse
import os

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, Distance, BarycentricMeanEcliptic
from sqlalchemy import text

from corgidb.ingest import gen_engine


POSITION_COLS = ["ra", "dec"]
KINEMATIC_COLS = ["sy_plx", "sy_pmra", "sy_pmdec", "st_radv"]
ALL_SAMPLED_COLS = POSITION_COLS + KINEMATIC_COLS

UNITS = {
    "ra": u.deg,
    "dec": u.deg,
    "sy_plx": u.mas,
    "sy_pmra": u.mas / u.yr,
    "sy_pmdec": u.mas / u.yr,
    "st_radv": u.km / u.s,
}

ERR_COLS = {
    "ra": ("raerr1", "raerr2"),
    "dec": ("decerr1", "decerr2"),
    "sy_plx": ("sy_plxerr1", "sy_plxerr2"),
    "sy_pmra": ("sy_pmraerr1", "sy_pmraerr2"),
    "sy_pmdec": ("sy_pmdecerr1", "sy_pmdecerr2"),
    "st_radv": ("st_radverr1", "st_radverr2"),
}

SQL_SELECT = """
    SELECT st_id, main_id, st_name,
           ra, raerr1, raerr2,
           `dec`, decerr1, decerr2,
           sy_plx, sy_plxerr1, sy_plxerr2,
           sy_pmra, sy_pmraerr1, sy_pmraerr2,
           sy_pmdec, sy_pmdecerr1, sy_pmdecerr2,
           st_radv, st_radverr1, st_radverr2
    FROM Stars
    WHERE (elat IS NULL OR elon IS NULL)
      AND ra IS NOT NULL AND `dec` IS NOT NULL
"""

SQL_UPDATE = """UPDATE Stars SET
    elat = :elat, elaterr1 = :elaterr1, elaterr2 = :elaterr2,
    elon = :elon, elonerr1 = :elonerr1, elonerr2 = :elonerr2
    WHERE st_id = :st_id"""


def two_piece_normal(rng, mu, sigma_upper, sigma_lower, n):
    """Draw n samples from a two-piece (split) normal distribution.

    The distribution has mode at mu, std sigma_upper for x > mu, and
    sigma_lower for x < mu.  Each side is drawn proportionally to its
    width so the mode (not the mean) is mu, which matches the convention
    used by err1/err2 columns reported as (upper, lower) one-sigma bounds.

    sigma_upper or sigma_lower of zero (or NaN) collapses that side to
    a delta function at mu; if both are zero, all samples equal mu.

    Args:
        rng (numpy.random.Generator): Source of randomness.
        mu (float): Mode of the distribution.
        sigma_upper (float): One-sigma half-width on the upper side.
        sigma_lower (float): One-sigma half-width on the lower side.
        n (int): Number of samples.

    Returns:
        numpy.ndarray: Length-n array of samples.
    """
    if not np.isfinite(sigma_upper) or sigma_upper < 0:
        sigma_upper = 0.0
    if not np.isfinite(sigma_lower) or sigma_lower < 0:
        sigma_lower = 0.0
    if sigma_upper == 0.0 and sigma_lower == 0.0:
        return np.full(n, mu, dtype=float)

    abs_z = np.abs(rng.standard_normal(n))
    total = sigma_upper + sigma_lower
    p_upper = sigma_upper / total
    upper = rng.random(n) < p_upper
    return mu + np.where(upper, abs_z * sigma_upper, -abs_z * sigma_lower)


def _safe_float(v):
    """Coerce a possibly-None/NaN scalar to float, returning NaN if invalid."""
    if v is None:
        return np.nan
    try:
        f = float(v)
    except (TypeError, ValueError):
        return np.nan
    return f


def sample_row(row, n, rng):
    """Build a length-n SkyCoord of MC samples for one star row.

    Position (ra, dec) is always sampled.  Kinematic args (sy_plx,
    sy_pmra, sy_pmdec, st_radv) are sampled when their central value is
    finite; if the central value is NaN, that argument is omitted from
    SkyCoord construction.  Per-axis NULL errors collapse to sigma=0 on
    that side, so passing n=1 with every error zero/NULL yields the
    central transform exactly.

    Args:
        row (pandas.Series): One DB row (must contain the columns listed
            in ALL_SAMPLED_COLS plus their *err1/*err2 partners).
        n (int): Number of MC samples.
        rng (numpy.random.Generator): Source of randomness.

    Returns:
        astropy.coordinates.SkyCoord:
            n samples already transformed to BarycentricMeanEcliptic.
    """
    samples = {}
    for col in ALL_SAMPLED_COLS:
        mu = _safe_float(row[col])
        if not np.isfinite(mu):
            continue
        e1, e2 = ERR_COLS[col]
        s_up = _safe_float(row.get(e1, np.nan))
        s_lo = _safe_float(row.get(e2, np.nan))
        s_up = 0.0 if not np.isfinite(s_up) else s_up
        s_lo = 0.0 if not np.isfinite(s_lo) else s_lo
        samples[col] = two_piece_normal(rng, mu, s_up, s_lo, n)

    kw = dict(frame="icrs", equinox="J2000", obstime="J2000")
    kw["ra"] = samples["ra"] * UNITS["ra"]
    kw["dec"] = samples["dec"] * UNITS["dec"]

    if "sy_plx" in samples:
        # Distance(parallax=...) rejects non-positive values; clip
        # negative draws to a tiny positive parallax so noisy samples
        # do not discard the rest of the kinematics for the row.
        plx = np.clip(samples["sy_plx"], 1e-6, None)
        kw["distance"] = Distance(parallax=plx * UNITS["sy_plx"])
    if "sy_pmra" in samples:
        kw["pm_ra_cosdec"] = samples["sy_pmra"] * UNITS["sy_pmra"]
    if "sy_pmdec" in samples:
        kw["pm_dec"] = samples["sy_pmdec"] * UNITS["sy_pmdec"]
    if "st_radv" in samples:
        kw["radial_velocity"] = samples["st_radv"] * UNITS["st_radv"]

    return SkyCoord(**kw).transform_to(BarycentricMeanEcliptic)


def _wrap_lon_for_percentile(lon_deg, central_lon_deg):
    """Re-center a wrapped longitude array around the central value.

    Longitudes naturally live on [0, 360); a star near the prime meridian
    will have samples on both sides of the wrap, and a naive percentile
    over those samples returns nonsense (~180 deg).  Shifting the array
    so the central value sits at 180 keeps the bulk of the distribution
    away from either edge before taking percentiles.
    """
    shifted = (lon_deg - central_lon_deg + 180.0) % 360.0
    return shifted, central_lon_deg


def _unshift_lon(percentile_value, central_lon_deg):
    return (percentile_value + central_lon_deg - 180.0) % 360.0


def compute_row(row, n, rng):
    """Compute (elat, elon, elaterr1, elaterr2, elonerr1, elonerr2) for one row.

    When all four ra/dec error columns are NULL the row falls back to
    the central transform only and the four error outputs are None.
    """
    pos_err_cols = ("raerr1", "raerr2", "decerr1", "decerr2")
    has_pos_err = all(np.isfinite(_safe_float(row[c])) for c in pos_err_cols)

    # Always evaluate the central transform so elat/elon are populated
    # even when error columns are missing.
    central = sample_row(row, 1, rng=np.random.default_rng(0))
    elat_c = float(central.lat[0].to_value(u.deg))
    elon_c = float(central.lon[0].to_value(u.deg))

    if not has_pos_err:
        return elat_c, elon_c, None, None, None, None

    ec = sample_row(row, n, rng)
    lat = ec.lat.to_value(u.deg)
    lon = ec.lon.to_value(u.deg)

    # Latitude is bounded on [-90, 90], percentiles are well-defined.
    p16, p50, p84 = np.percentile(lat, [16, 50, 84])
    elat = float(p50)
    elaterr1 = float(p84 - p50)
    elaterr2 = float(p50 - p16)

    # Longitude is circular; re-center around the central value before
    # taking percentiles so the wrap at 0/360 deg does not corrupt them.
    lon_shift, _ = _wrap_lon_for_percentile(lon, elon_c)
    p16, p50, p84 = np.percentile(lon_shift, [16, 50, 84])
    elon = float(_unshift_lon(p50, elon_c))
    elonerr1 = float(p84 - p50)
    elonerr2 = float(p50 - p16)

    return elat, elon, elaterr1, elaterr2, elonerr1, elonerr2


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--username", default=os.getenv("PLANDB_USER"))
    p.add_argument("--db", default=os.getenv("PLANDB_NAME", "plandb_scratch"))
    p.add_argument("--server", default=os.getenv("PLANDB_SERVER", "127.0.0.1"))
    p.add_argument(
        "--n-samples",
        type=int,
        default=10000,
        dest="n_samples",
        help="Monte Carlo samples per row (default: 10000).",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
        help="PRNG seed for reproducibility (default: 42).",
    )
    p.add_argument(
        "--write",
        action="store_true",
        help="Commit results to the DB (default: dry run).",
    )
    return p.parse_args()


def main(args):
    if not args.username:
        raise SystemExit("Missing DB username; pass --username or set PLANDB_USER.")

    print(f"Using database: {args.db} on {args.server}")
    print("Mode:", "WRITE" if args.write else "DRY RUN")
    if not args.write:
        print("Re-run with --write to commit changes.")

    engine = gen_engine(args.username, db=args.db, server=args.server)

    with engine.connect() as conn:
        df = pd.read_sql(text(SQL_SELECT), conn)

    print(f"\nFound {len(df)} row(s) with missing elat/elon (and ra/dec present).")
    if df.empty:
        return

    rng = np.random.default_rng(args.seed)

    rows_out = []
    n_with_errors = 0
    n_central_only = 0
    for _, row in df.iterrows():
        elat, elon, e1lat, e2lat, e1lon, e2lon = compute_row(row, args.n_samples, rng)
        if e1lat is None:
            n_central_only += 1
        else:
            n_with_errors += 1
        rows_out.append(
            {
                "st_id": int(row.st_id),
                "elat": elat,
                "elaterr1": e1lat,
                "elaterr2": e2lat,
                "elon": elon,
                "elonerr1": e1lon,
                "elonerr2": e2lon,
            }
        )

    print(f"  With MC uncertainties : {n_with_errors} row(s)")
    print(f"  Central value only    : {n_central_only} row(s)")

    out_df = pd.DataFrame(rows_out)

    if not args.write:
        with pd.option_context("display.width", 160, "display.max_columns", None):
            print("\nSample of computed values (first 20):")
            print(out_df.head(20).to_string(index=False))
        print("\nDry run; no rows written.")
        return

    with engine.begin() as conn:
        for r in rows_out:
            conn.execute(text(SQL_UPDATE), r)
    print(f"\nWrote {len(rows_out)} row(s) to Stars.")


if __name__ == "__main__":
    main(parse_args())
