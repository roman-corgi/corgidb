import numpy as np
from sqlalchemy import text
from sqlalchemy.exc import SQLAlchemyError
from astroquery.vizier import Vizier
from astroquery.exceptions import TableParseError 
from astropy.coordinates import SkyCoord
from astropy import units as u
from sqlalchemy import create_engine
import os
import sys
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # .../corgidbcopy
CORGI_DIR = os.path.join(ROOT_DIR, "corgidb")
sys.path.append(CORGI_DIR)
from ingest import gen_engine   

# connect to db  (update with db=plandb when done testing)
engine = gen_engine(username="plandb_user", db="plandb_scratch", server="127.0.0.1")

#change last line of code when done testing

Vizier.ROW_LIMIT = 50

def _find_jsdc_jmdc():
    jsdc = list(Vizier.find_catalogs("JSDC").keys())
    jmdc = list(Vizier.find_catalogs("JMDC").keys())
    return jsdc, jmdc

_JSDC_IDS, _JMDC_IDS = _find_jsdc_jmdc()

def _extract_theta(table):
    if len(table) == 0:
        return None

    cols = table.colnames

    #JSDC: limb-darkened diameter
    if 'AngDiam' in cols:
        theta = table['AngDiam'][0]
        if theta is not None and not np.isnan(theta):
            if 'e_AngDiam' in cols:
                error = table['e_AngDiam'][0]
            else:
                error = 0.1 * theta
            return float(theta), float(error)

    #JMDC: prefer limb-darkened (LDD), else uniform-disk (UDD)
    for main_col, err_col in [('LDD', 'e_LDD'), ('UDD', 'e_UDD')]:
        if main_col in cols:
            theta = table[main_col][0]
            if theta is None or np.isnan(theta):
                continue
            if err_col in cols:
                error = table[err_col][0]
            else:
                error = 0.1 * theta
            return float(theta), float(error)

    #generic column names for ang diameter
    possible_cols = ['theta', 'diameter', 'ang_diam', 'D_']
    for col in cols:
        if col in possible_cols:
            theta = table[col][0]
            if theta is not None and not np.isnan(theta):
                # try two common error conventions: col_err, e_col
                for err_col in (col + '_err', 'e_' + col):
                    if err_col in cols:
                        error = table[err_col][0]
                        break
                else:
                    error = 0.1 * theta
                return float(theta), float(error)

    return None

#find stellar diameter from JDSC or JMDC
def _query_vizier_ra_dec(ra_deg, dec_deg, radius_arcsec=60.0):
    if ra_deg is None or dec_deg is None:
        return None
    coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame="icrs")

    #JSDC first
    for cat in _JSDC_IDS:
        try:
            res = Vizier.query_region(coord, radius=radius_arcsec*u.arcsec, catalog=cat)
            if res and len(res) > 0 and len(res[0]) > 0:
                got = _extract_theta(res[0])
                if got:
                    t, e = got
                    return t, e, f"JSDC:{cat}"
        except (TableParseError, Exception):
            pass

    #JMDC
    for cat in _JMDC_IDS:
        try:
            res = Vizier.query_region(coord, radius=radius_arcsec*u.arcsec, catalog=cat)
            if res and len(res) > 0 and len(res[0]) > 0:
                got = _extract_theta(res[0])
                if got:
                    t, e = got
                    return t, e, f"JMDC:{cat}"
        except (TableParseError, Exception):
            pass

    return None


#exosims stellar diameter calc
def calculate_stellar_diameter(BV, Vmag):
    """Approximate stellar angular diameter in mas using BV color and Vmag.

    Computed according to the model from [Boyajian2014]_. Returns diameter only.
    """
    if BV is None or Vmag is None:
        return None

    # B-V polynomial fit coefficients:
    coeffs = [0.49612, 1.11136, -1.18694, 0.91974, -0.19526]

    # evaluate eq. 2 using B-V color
    logth_zero = 0.0
    for j, ai in enumerate(coeffs):
        logth_zero += ai * BV ** j

    # now invert eq. 1 to get the actual diameter in mas
    theta = 10 ** (logth_zero - 0.2 * Vmag)
    return float(theta)



def main(dry_run=True, limit=10):
    # identifies stars missing diameters 
    with engine.begin() as con:
        result = con.execute(text("SHOW TABLES;"))
        print("Tables in DB:")
        for row in result:
            print(row)

        result = con.execute(text("SHOW COLUMNS FROM Stars;"))
        print("\nColumns in Stars table:")
        columns = [row['Field'] for row in result]
        print(columns)

        query = text(f"""
            SELECT
                st_id      AS id,
                st_name    AS name,
                CAST(ra   AS DECIMAL(20,10))  AS ra_deg,
                CAST(`dec` AS DECIMAL(20,10)) AS dec_deg,
                (sy_bmag - sy_vmag) AS BV,
                sy_vmag    AS Vmag
            FROM Stars
            WHERE diameter IS NULL
            LIMIT {limit}
        """)
        result = con.execute(query)
        missing_stars = result.fetchall()
        print(f"\nFound {len(missing_stars)} Stars with missing diameters.")

    stars_to_update = [
        dict(zip(['id', 'name', 'ra_deg', 'dec_deg', 'BV', 'Vmag'], row))
        for row in missing_stars
    ]
    print(f"Sample star: {stars_to_update[0] if stars_to_update else 'None'}")

    #Compute diameters (Vizier or fallback calc)
    updates = []
    for star in stars_to_update:
        star_id = star['id']
        star_name = star['name']
        ra = star['ra_deg']
        dec = star['dec_deg']
        bv = star['BV']
        vmag = star['Vmag']

        print(f"\nProcessing {star_name} (ID: {star_id})")
         # ensure ra/dec are usable numbers
        try:
            ra = float(ra) if ra is not None else None
            dec = float(dec) if dec is not None else None
        except (TypeError, ValueError):
            print(f"  Skipping {star_name}: bad RA/Dec ({ra}, {dec})")
            continue

        vizier_result = _query_vizier_ra_dec(ra, dec)
        if vizier_result:
            theta, error, source = vizier_result
            print(f"  Got from {source}: {theta:.3f} ± {error:.3f} mas")
        else:
            theta = calculate_stellar_diameter(bv, vmag)
            if theta is not None:
                source = "Calculated (BV/Vmag)"
                print(f"  Calculated: {theta:.3f} mas")
            else:
                print("  Failed: No RA/Dec or BV/Vmag available")
                continue

        updates.append({'id': star_id, 'diameter': theta})

    #dry run vs actual update
    if dry_run:
        print("\nDRY RUN – would update the following rows:")
        for u in updates:
            print(f"  ID {u['id']} -> diameter = {u['diameter']:.3f} mas")
        return

    success_count = 0
    with engine.begin() as con:
        for update in updates:
            try:
                query = text("""
                    UPDATE Stars 
                    SET diameter = :diameter 
                    WHERE st_id = :id
                """)
                con.execute(query, {
                    'diameter': update['diameter'],
                    'id': update['id']
                })
                success_count += 1
                print(f"Updated ID {update['id']}: {update['diameter']:.3f} mas")
            except SQLAlchemyError as e:
                print(f"Error updating ID {update['id']}: {e}")

    print(f"\ncomplete: Updated {success_count}/{len(updates)} Stars.")


if __name__ == "__main__":
    main(dry_run=True, limit=50) 
    #first test only prints (current # of stars updated at a time: limit is currently 50)
    #change to False to make changes (when done testing)
    #change limit to update more stars at a time