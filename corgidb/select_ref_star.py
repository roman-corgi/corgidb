"""Module For Refrence Star selection and supporting funcitons"""
import numpy as np
import astropy.units as u
import astropy.time as t
import astropy.coordinates as c
import pandas as pd
import sqlalchemy as sql
from roman_pointing import roman_pointing as rp

# Wrapper function to return a dataframe rather than an sql object and clean up the sqlalchemy objects inside
def select_query_db(conn: sql.engine.base.Connection, stmt: sql.Select) -> pd.DataFrame:
    """Take db connection and select statment to return data from the db

    Args:
        conn (sqlalchemy.engine.base.Connection): sqlalchemy connection object
        stmt (sqlalchemy.Select): sqlalchemy select statement for desired data

    Returns:
        data (pandas.DataFrame): dataframe containing the desired data
    """
    #change return into a dataframe
    raw_data = pd.read_sql(stmt, conn)
    #select only data thats not sqlalchemy objects and return them
    good_data = ~raw_data.columns.str.startswith('_sa_')
    data = raw_data.loc[:, good_data]
    return data

def check_pointing(tar: pd.DataFrame, obs_start: t.Time, obs_duration: t.TimeDelta, slices=100) -> tuple[bool, list, list]:
    """Check that the observation window defined for a stars data does not violate any constraints
    
    Args:
        tar (pandas.DataFrame): Target star data defined as a single row DataFrame
        obs_start (astropy.time.Time): Start time of the observation
        obs_duration (astropy.time.Time): Duration of the observation
        slices (int): The number of time slices to calculate angles for
    
    Returns:
        result (tuple[bool, list, list]): returns a boolean true/false if the observation is valid, and a list of pointings over the window for sun and pitch angles
    """
    #Calculate times array to calculate each angle at.
    times = obs_start + obs_duration * np.linspace(0, 1, slices)
    #Create Skycoord object for the target at the start time 
    tar_cords = c.SkyCoord(
    tar.loc[0, "ra"] * u.degree,
    tar.loc[0, "dec"] * u.degree,
    unit=(u.degree, u.degree),
    frame="icrs",
    distance=tar.loc[0,"sy_dist"] * u.parsec, 
    pm_ra_cosdec=tar.loc[0,"sy_pmra"] * u.milliarcsecond / u.year, 
    pm_dec=tar.loc[0,"sy_pmdec"] * u.milliarcsecond / u.year,
    radial_velocity=tar.loc[0, "st_radv"] * u.km / u.second,
    equinox="J2000",
    obstime="J2000",
    ).transform_to(c.BarycentricMeanEcliptic)
    #calculate angles of interest over the observation window 
    sun_ang_targ, yaw_targ, pitch_targ, B_C_I_targ = rp.calcRomanAngles(
    tar_cords, times, rp.getL2Positions(times)
    )
    #Convert angles to degrees
    sun_ang_d_targ = sun_ang_targ.to(u.degree)
    pitch_d_targ = pitch_targ.to(u.degree)
    #Check the five degree sun angle constraint and set validity
    if any((item < 54 * u.degree) or (item > 126 * u.degree) for item in sun_ang_d_targ):
        valid = False
    else:
        valid = True
    result = valid, sun_ang_d_targ, pitch_d_targ
    return result

def select_ref_star(st_name: str, obs_start: t.Time, obs_duration: t.TimeDelta, engine: sql.engine.base.Engine) -> str:
    """Select refrence star given target and observation parameters

    Args: 
        st_name (str): Target star name in db
        obs_start (astropy.time.Time): Start time of observation window
        obs_duration (astropy.time.Time): Duration of the observation window
        engine (sql.engine.base.Engine): Sqlalchemy engine object that is connected to plandb

    Returns:
        str: Reference star name in db
    """
    # Type checking for inputs and other generic error handling

    # Connect to the DB and get tables
    metadata = sql.MetaData()
    stars_table = sql.Table('Stars', metadata, autoload_with=engine)
    # pass this connection to the query method to be used and then spun up and cleaned up in the calling method 
    conn = engine.connect()
    # query for a Star with the correct st_name entry
    stmt = sql.select(stars_table).where(stars_table.c.st_name == st_name)
    # Test values/statements for making sure that little bits of the code work
    sci_target = select_query_db(conn, stmt)
    tar_val, tar_sun_angs, tar_pitch_angs = check_pointing(sci_target, obs_start, obs_duration)
    if tar_val is False:
        ref_star = "Observation window violates Solar angle Constraint"
        print(ref_star)
        # logic to handle a observation that is not valid
    else:
        ref_star = ""
        ref_grade = ['A','B','C']
        i = 0
        while ref_star is "" and i < 3:
            cur_ref_grade = ref_grade[i]
            ref_stmt= sql.select(stars_table).where(stars_table.c.st_psfgrade == cur_ref_grade)
            ref_stars_data = select_query_db(conn, ref_stmt)
            for record in ref_stars_data.to_dict(orient="records"):
                cur_tar_star = pd.DataFrame.from_dict(record)
                ref_val, ref_sun_angs, ref_pitch_angs = check_pointing(cur_tar_star, obs_start, obs_duration)
                if ref_val:
                    del_pitch = np.array(tar_pitch_angs) - np.array(ref_pitch_angs)
                    max_pitch = np.abs(np.max(del_pitch))
                    if max_pitch < 5:
                        ref_star = record.loc[0, "st_name"]
                        i = 3
            i = i + 1
        if ref_star == "":
            ref_star = f"No refrence star of class A, B, or C was found for {st_name} between {obs_start} and {obs_start+obs_duration}"
        # Logic to start selecting reference stars
    return ref_star