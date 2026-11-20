"""SIMBAD and VizieR mirror selection.

All SIMBAD access in corgidb goes through :func:`get_simbad`, which tries each
mirror in :data:`SIMBAD_SERVERS` in order and returns ``None`` (after printing a
message) when none of them respond. :data:`VIZIER_SERVERS` lists the VizieR
mirrors in the same priority order: CDS, then the Harvard CfA mirror, then the
South Africa IDIA mirror.
"""

from typing import Optional

import requests
from astroquery.simbad import SimbadClass
from astroquery.simbad import conf as simbad_conf
from pyvo.dal.exceptions import DALFormatError, DALServiceError

SIMBAD_SERVERS = ("simbad.cds.unistra.fr", "simbad.harvard.edu", "simbad.idia.ac.za")
VIZIER_SERVERS = (
    "vizier.cds.unistra.fr",
    "vizier.cfa.harvard.edu",
    "vizier.idia.ac.za",
)
SIMBAD_ERRORS = (DALServiceError, DALFormatError, requests.RequestException)


def get_simbad(*votable_fields: str) -> Optional[SimbadClass]:
    """Return a Simbad instance bound to the first reachable SIMBAD mirror.

    Mirrors are tried in the order given by :data:`SIMBAD_SERVERS`. Each
    candidate is probed with a TAP capabilities request before the requested
    VOTable fields are added.

    Args:
        *votable_fields (str):
            VOTable fields to add to the returned instance (e.g. ``"pmra"``).

    Returns:
        Optional[astroquery.simbad.SimbadClass]:
            Configured Simbad instance, or None if no mirror is available.

    Examples:
        >>> simbad = get_simbad("pmra", "pmdec")
        >>> if simbad is not None:
        ...     res = simbad.query_object("HD 95128")

    """

    for server in SIMBAD_SERVERS:
        # astroquery only accepts servers in this list; it does not know about
        # the IDIA mirror.
        if server not in simbad_conf.servers_list:
            simbad_conf.servers_list.append(server)
        simbad = SimbadClass()
        simbad.server = server
        try:
            simbad.hardlimit
            simbad.add_votable_fields(*votable_fields)
        except SIMBAD_ERRORS:
            print(f"SIMBAD mirror {server} unreachable; trying next.")
            continue
        return simbad

    print("No SIMBAD mirror is available.")
    return None
