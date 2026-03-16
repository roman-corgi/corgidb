import requests
import pandas

url = "https://corgidb.sioslab.com/fetch_refs.php"
response = requests.get(url, headers={"User-Agent": "XY"})

if response.status_code == 200:
    data = response.json()
else:
    print("Could not retrieve data")

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
