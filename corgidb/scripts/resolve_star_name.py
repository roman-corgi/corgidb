from typing import Dict, Iterable, Union

import requests


def resolve_star_name(
    st_names: Union[str, Iterable[str]],
    url: str = "https://corgidb.sioslab.com/resolve_star_name.php",
) -> Dict[str, str]:
    """Resolve one or more star names/aliases to their main_id.

    Resolution happens via the StarAliases table. A single request is made
    regardless of how many names are passed in.

    Args:
        st_names (Union[str, Iterable[str]]):
            Star name/alias, or names/aliases, to resolve (e.g. ``"47 UMa"``
            or ``["47 UMa", "HD 95128"]``).
        url (str):
            URL of the resolve_star_name.php endpoint.

    Returns:
        Dict[str, str]:
            Mapping of input name to resolved main_id. Names with no match
            in StarAliases are omitted.

    """

    if isinstance(st_names, str):
        st_names = [st_names]
    else:
        st_names = list(st_names)

    response = requests.get(
        url,
        headers={"User-Agent": "corgidb_agent"},
        params={"st_names": ",".join(st_names)},
    )

    response.raise_for_status()

    data = response.json()

    return {alias: main_id for alias, main_id in data}
