import requests
import pandas
import numpy as np

url = "https://corgidb.sioslab.com/fetch_star.php"
response = requests.get(
    url, headers={"User-Agent": "corgidb_agent"}, params={"st_name": "47 UMa"}
)

assert response.status_code == 200, "Query failed."

data = response.json()
data = np.vstack(data).transpose()

colnames = [
    "st_name",
    "main_id",
    "ra",
    "dec",
    "spectype",
    "sy_vmag",
    "sy_imag",
    "sy_dist",
    "sy_plx",
    "sy_pmra",
    "sy_pmdec",
    "st_radv",
]

out = {}
for colname, col in zip(colnames, data):
    out[colname] = col

out = pandas.DataFrame(out)
