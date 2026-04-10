import numpy as np
import pandas
import requests


def query_star(
    st_name: str,
    url: str = "https://corgidb.sioslab.com/fetch_star.php",
) -> pandas.DataFrame:
    """Query a star by name or alias.

    The server resolves aliases via the StarAliases table before
    querying the Stars table, so both primary names and aliases work.

    Args:
        st_name (str):
            Star name or alias to look up (e.g. ``"47 UMa"``,
            ``"HD 95128"``).
        url (str):
            URL of the fetch_star.php endpoint.

    Returns:
        pandas.DataFrame:
            Query results with columns: st_name, main_id, ra, dec,
            spectype, sy_vmag, sy_imag, sy_dist, sy_plx, sy_pmra,
            sy_pmdec, st_radv. Empty DataFrame if no match found.

    """

    response = requests.get(
        url, headers={"User-Agent": "corgidb_agent"}, params={"st_name": st_name}
    )

    response.raise_for_status()

    data = response.json()

    colnames: list[str] = [
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

    if len(data) == 0:
        return pandas.DataFrame(columns=colnames)

    data = np.vstack(data).transpose()

    out: dict[str, np.ndarray] = {}
    for colname, col in zip(colnames, data):
        out[colname] = col

    return pandas.DataFrame(out)
