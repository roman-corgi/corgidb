import numpy as np
import astropy.units as u
import astropy.time as t
import sqlalchemy as sql
from corgidb import select_ref_star as cdb



st_name = "47 Uma"
eng = cdb.ingest.gen_engine('plandb_user', 'plandb_scratch')
metadata = sql.MetaData()
stars_table = sql.Table('Stars', metadata, autoload_with=eng)
conn = eng.connect()
stmt = sql.select(stars_table).where(stars_table.c.st_name == st_name)
data = cdb.select_ref_star.select_query_db(conn, stmt)
print(data["dec"])
t_str = ["2027-01-01T00:00:00.0"]
o_s = t.Time(t_str, format="isot", scale="utc")
o_d = 365 * u.day
val, sun_ang, pitch_ang = check_pointing(data, o_s, o_d)
t = np.linspace(0,365,100)
print(val)
print(sun_ang)
plt.plot(t, sun_ang)
plt.show()
print(pitch_ang)