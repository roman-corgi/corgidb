import requests

url = "https://corgidb.sioslab.com/fetch_refs.php"
response = requests.get(url, headers={"User-Agent": "XY"})

if response.status_code == 200:
    data = response.json()
    for d in data:
        print(d)
