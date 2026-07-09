"""Example: Filter and sort reference stars by observable properties.

Shows common filtering patterns:
  - Stars brighter than a V-band magnitude limit
  - Stars with a measured uniform disk diameter
  - Stars sorted by PSF grade for a specific bandpass
"""

import pandas as pd
from corgidb import CorgiQuery

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 140)

cq = CorgiQuery()
refs = cq.query_refstars()

# Cast numeric columns (API returns strings)
numeric_cols = [
    "sy_vmag", "sy_imag", "sy_dist",
    "st_uddv", "st_uddi",
    "st_psfgrade_nfb1_high", "st_psfgrade_nfb1_med",
    "st_psfgrade_specb3_high", "st_psfgrade_specb3_med",
    "st_psfgrade_wfb4_high", "st_psfgrade_wfb4_med",
]
for col in numeric_cols:
    refs[col] = pd.to_numeric(refs[col], errors="coerce")

# --- Filter 1: V-band magnitude brighter than 7 ---
bright = refs[refs["sy_vmag"] < 7].copy()
print(f"Stars with V < 7: {len(bright)}")
print(bright[["st_name", "spectype", "sy_vmag", "sy_dist"]].to_string(index=False))
print()

# --- Filter 2: stars with a measured V-band uniform disk diameter ---
has_uddv = refs[refs["st_uddv"].notna() & (refs["st_uddv"] > 0)].copy()
print(f"Stars with measured ud diameter (V-band): {len(has_uddv)}")
print(has_uddv[["st_name", "spectype", "sy_vmag", "st_uddv"]].to_string(index=False))
print()

# --- Filter 3: best PSF references for NFB1 high-contrast ---
grade_col = "st_psfgrade_nfb1_high"
top_psf = refs[refs[grade_col].notna()].sort_values(grade_col, ascending=False).head(10)
print(f"Top 10 reference stars by {grade_col}:")
print(top_psf[["st_name", "sy_vmag", grade_col]].to_string(index=False))
