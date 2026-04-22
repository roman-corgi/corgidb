import pandas

from corgidb.scripts.query_refstars import query_refstars
from corgidb.scripts.query_star import query_star

_DEFAULT_BASE_URL = "https://corgidb.sioslab.com"


class CorgiQuery:
    """Client for querying the corgidb database via its HTTP API.

    Wraps the query functions in ``corgidb.scripts`` to provide a
    session-like interface with a configurable base URL.

    Star name queries automatically resolve aliases through the
    StarAliases table on the server, so you may query by any
    registered alias (e.g. ``"HD 95128"`` instead of ``"47 UMa"``).

    Args:
        base_url (str):
            Root URL of the corgidb frontend. Defaults to the production
            deployment at ``https://corgidb.sioslab.com``.

    """

    def __init__(self, base_url: str = _DEFAULT_BASE_URL) -> None:
        self.base_url = base_url.rstrip("/")

    def query_star(self, st_name: str) -> pandas.DataFrame:
        """Query a star by name or alias.

        Args:
            st_name (str):
                Star name or alias to look up.

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
                Empty DataFrame if no match found.

        """

        return query_star(st_name, url=f"{self.base_url}/fetch_star.php")

    def query_refstars(self) -> pandas.DataFrame:
        """Query all reference stars.

    Returns:
        pandas.DataFrame:
            All reference stars with columns:
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
                st_radv,
                st_psfgrade_nfb1_high,
                st_psfgrade_nfb1_med,
                st_psfgrade_specb3_high,
                st_psfgrade_specb3_med,
                st_psfgrade_wfb4_high,
                st_psfgrade_wfb4_med,
                st_uddv,
                st_uddi,
                st_uddmeas,
                st_lddmeas.
            Empty DataFrame if no results.

        """

        return query_refstars(url=f"{self.base_url}/fetch_refs.php")
