import numpy as np
import pandas
import requests


def query_refstars(
    url: str = "https://corgidb.sioslab.com/fetch_refs.php",
) -> pandas.DataFrame:
    """Query all reference stars.

    Args:
        url (str):
            URL of the fetch_refs.php endpoint.

    Returns:
        pandas.DataFrame:
            All reference stars with columns: st_name, main_id, ra, dec,
            spectype, sy_vmag, sy_imag, sy_dist, sy_plx, sy_pmra, sy_pmdec,
            st_radv, st_psfgrade_nfb1_high, st_psfgrade_nfb1_med,
            st_psfgrade_specb3_high, st_psfgrade_specb3_med,
            st_psfgrade_wfb4_high, st_psfgrade_wfb4_med, st_uddv, st_uddi,
            st_uddmeas, st_lddmeas. Empty DataFrame if no results.

    """

    response = requests.get(url, headers={"User-Agent": "corgidb_agent"})

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
        "st_psfgrade_nfb1_high",
        "st_psfgrade_nfb1_med",
        "st_psfgrade_specb3_high",
        "st_psfgrade_specb3_med",
        "st_psfgrade_wfb4_high",
        "st_psfgrade_wfb4_med",
        "st_uddv",
        "st_uddi",
        "st_uddmeas",
        "st_lddmeas",
    ]

    if len(data) == 0:
        return pandas.DataFrame(columns=colnames)

    data = np.vstack(data).transpose()

    out: dict[str, np.ndarray] = {}
    for colname, col in zip(colnames, data):
        out[colname] = col

    return pandas.DataFrame(out)
