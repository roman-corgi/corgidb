"""Example: Query the full reference star catalog.

Retrieves all reference stars and shows how to inspect PSF grades
and uniform disk diameter columns.
"""

import pandas as pd
from corgidb import CorgiQuery

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 120)

cq = CorgiQuery()
refs = cq.query_refstars()

print(f"Total reference stars: {len(refs)}")
print()

# Basic stellar info
print("--- Stellar parameters (first 10 rows) ---")
cols_basic = ["st_name", "main_id", "spectype", "sy_vmag", "sy_imag", "sy_dist"]
print(refs[cols_basic].head(10).to_string(index=False))
print()

# PSF grade columns
psf_cols = [c for c in refs.columns if c.startswith("st_psfgrade")]
print("--- PSF grades (first 10 rows) ---")
print(refs[["st_name"] + psf_cols].head(10).to_string(index=False))
print()

# Uniform disk diameter columns
ud_cols = ["st_uddv", "st_uddi", "st_uddmeas", "st_lddmeas"]
print("--- Disk diameters (first 10 rows) ---")
print(refs[["st_name"] + ud_cols].head(10).to_string(index=False))
