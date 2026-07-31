from typing import List

import numpy as np
import pandas
import requests

from corgidb.scripts.columns import STAR_COLUMNS


def query_stars(
    st_names: List[str],
    url: str = "https://corgidb.sioslab.com/fetch_stars.php",
) -> pandas.DataFrame:
    """Query multiple stars by name or alias.

    The server resolves aliases via the StarAliases table before
    querying the Stars table, so both primary names and aliases work.
    Names that don't resolve to a known star are silently omitted from
    the result.

    Args:
        st_names (List[str]):
            Star names or aliases to look up (e.g. ``["47 UMa", "HD 95128"]``).
        url (str):
            URL of the fetch_stars.php endpoint.

    Returns:
        pandas.DataFrame:
            Query results with columns:
                st_name,
                main_id,
                ra,
                dec,
                spectype,
                sy_vmag,
                sy_imag,
                sy_dist,
                sy_plx,
                sy_pmra,
                sy_pmdec,
                st_radv.
            Empty DataFrame if no matches found.

    """

    response = requests.get(
        url,
        headers={"User-Agent": "corgidb_agent"},
        params={"st_names": ",".join(st_names)},
    )

    response.raise_for_status()

    data = response.json()

    if len(data) == 0:
        return pandas.DataFrame(columns=STAR_COLUMNS)

    data = np.vstack(data).transpose()

    out: dict[str, np.ndarray] = {}
    for colname, col in zip(STAR_COLUMNS, data):
        out[colname] = col

    return pandas.DataFrame(out)
