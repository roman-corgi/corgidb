import os
import urllib.parse
from typing import Iterable, List, Union

import numpy as np
import pandas
from astroquery.simbad import Simbad

from corgidb.scripts.columns import STAR_COLUMNS
from corgidb.scripts.query_refstars import query_refstars
from corgidb.scripts.query_star import query_star
from corgidb.scripts.query_stars import query_stars
from corgidb.scripts.resolve_star_name import resolve_star_name
from corgidb.utils import get_cache_dir

_DEFAULT_BASE_URL = "https://corgidb.sioslab.com"


def _atomic_write_csv(df: pandas.DataFrame, path: str) -> None:
    """Write a DataFrame to CSV atomically.

    Args:
        df (pandas.DataFrame):
            DataFrame to write.
        path (str):
            Destination path.

    """

    tmp_path = f"{path}.tmp"
    df.to_csv(tmp_path, index=False)
    os.replace(tmp_path, path)


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

    Attributes:
        starcache (pandas.DataFrame):
            In-memory copy of the on-disk star cache used by
            :meth:`cache_stars`, keyed by ``main_id``. Loaded from disk on
            construction; empty if no cache file exists yet for this
            instance's ``base_url``.

    """

    def __init__(self, base_url: str = _DEFAULT_BASE_URL) -> None:
        self.base_url = base_url.rstrip("/")

        netloc = urllib.parse.urlparse(self.base_url).netloc
        self._starcache_path = os.path.join(get_cache_dir(), f"starcache_{netloc}.csv")
        if os.path.isfile(self._starcache_path):
            try:
                self.starcache = pandas.read_csv(
                    self._starcache_path, dtype={"main_id": str}
                )
            except (pandas.errors.EmptyDataError, pandas.errors.ParserError):
                self.starcache = pandas.DataFrame(columns=STAR_COLUMNS)
        else:
            self.starcache = pandas.DataFrame(columns=STAR_COLUMNS)

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

    def query_stars(self, st_names: List[str]) -> pandas.DataFrame:
        """Query multiple stars by name or alias.

        Args:
            st_names (List[str]):
                Star names or aliases to look up.

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

        return query_stars(st_names, url=f"{self.base_url}/fetch_stars.php")

    def _merge_into_cache(self, fetched: pandas.DataFrame) -> None:
        """Merge freshly-fetched rows into starcache and persist to disk.

        Any existing cache rows sharing a main_id with ``fetched`` are
        replaced, so ``starcache`` never holds duplicate main_id rows.

        Args:
            fetched (pandas.DataFrame):
                Newly-fetched rows to merge in. No-op if empty.

        """

        if fetched.empty:
            return

        self.starcache = pandas.concat(
            [
                self.starcache[~self.starcache["main_id"].isin(fetched["main_id"])],
                fetched,
            ],
            ignore_index=True,
        )
        _atomic_write_csv(self.starcache, self._starcache_path)

    def query_simbad(self, st_name: str) -> pandas.DataFrame:
        """Query SIMBAD directly for a star not found in corgidb.

        Used by :meth:`cache_stars` as a fallback for names that don't
        resolve via corgidb's StarAliases table. Only position, proper
        motion, parallax, and radial velocity are available from this
        query; spectype, sy_vmag, sy_imag, and sy_dist are set to null.

        Args:
            st_name (str):
                Star name to look up in SIMBAD.

        Returns:
            pandas.DataFrame:
                Single-row result with the same columns as query_stars.
                Empty DataFrame if SIMBAD has no match or the query fails.

        """

        simbad = Simbad()
        simbad.add_votable_fields("pmra", "pmdec", "plx_value", "rvz_radvel")

        try:
            res = simbad.query_object(st_name)
        except Exception:
            res = None

        if res is None or len(res) == 0:
            return pandas.DataFrame(columns=STAR_COLUMNS)

        for col in ["pmra", "pmdec", "rvz_radvel"]:
            res[col] = res[col].filled(0)

        row = res[0]
        plx = row["plx_value"]
        plx = None if np.ma.is_masked(plx) else float(plx)

        return pandas.DataFrame(
            [
                {
                    "st_name": st_name,
                    "main_id": str(row["main_id"]),
                    "ra": float(row["ra"]),
                    "dec": float(row["dec"]),
                    "spectype": None,
                    "sy_vmag": None,
                    "sy_imag": None,
                    "sy_dist": None,
                    "sy_plx": plx,
                    "sy_pmra": float(row["pmra"]),
                    "sy_pmdec": float(row["pmdec"]),
                    "st_radv": float(row["rvz_radvel"]),
                }
            ],
            columns=STAR_COLUMNS,
        )

    def cache_stars(
        self, st_names: Union[str, Iterable[str]], force_new: bool = False
    ) -> pandas.DataFrame:
        """Query stars by name or alias, using a local on-disk cache.

        Each name is resolved to its main_id via the StarAliases table.
        Names already present in the cache are served from disk; names not
        yet cached are fetched via :meth:`query_stars` and added to both the
        in-memory ``starcache`` and the on-disk cache file. Names that don't
        resolve via corgidb's StarAliases table fall back to a direct
        SIMBAD lookup (see :meth:`query_simbad`), which leaves spectype,
        sy_vmag, sy_imag, and sy_dist null. Names not found in either source
        are silently omitted, and duplicate aliases for the same star
        collapse to a single row.

        Args:
            st_names (Union[str, Iterable[str]]):
                Star name/alias, or names/aliases, to look up.
            force_new (bool):
                If True, ignore the cache for the requested names, always
                re-fetch fresh data, and overwrite any matching entries in
                the on-disk cache. Defaults to False.

        Returns:
            pandas.DataFrame:
                Query results, in first-occurrence order of the resolved
                input names, with columns:
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
                Empty DataFrame if no names resolve.

        """

        if isinstance(st_names, str):
            st_names = [st_names]
        else:
            st_names = list(st_names)

        resolved = resolve_star_name(
            st_names, url=f"{self.base_url}/resolve_star_name.php"
        )

        main_id_to_name: dict[str, str] = {}
        for name, main_id in resolved.items():
            main_id_to_name.setdefault(main_id, name)

        corgidb_ids: List[str] = list(dict.fromkeys(resolved.values()))

        cached_ids = set(self.starcache["main_id"])
        if force_new:
            to_fetch_ids = corgidb_ids
        else:
            to_fetch_ids = [i for i in corgidb_ids if i not in cached_ids]

        if to_fetch_ids:
            fetch_names = [main_id_to_name[i] for i in to_fetch_ids]
            self._merge_into_cache(self.query_stars(fetch_names))

        # SIMBAD fallback for names corgidb has never heard of.
        cached_ids = set(self.starcache["main_id"])
        simbad_ids_for_name: dict[str, str] = {}
        new_simbad_rows = []
        for name in st_names:
            if name in resolved:
                continue
            if not force_new and name in cached_ids:
                # already cached from an earlier SIMBAD lookup under this
                # exact name
                simbad_ids_for_name[name] = name
                continue
            simbad_row = self.query_simbad(name)
            if simbad_row.empty:
                continue
            simbad_ids_for_name[name] = simbad_row.iloc[0]["main_id"]
            new_simbad_rows.append(simbad_row)

        if new_simbad_rows:
            self._merge_into_cache(pandas.concat(new_simbad_rows, ignore_index=True))

        ordered_ids: List[str] = []
        for name in st_names:
            if name in resolved:
                mid = resolved[name]
            elif name in simbad_ids_for_name:
                mid = simbad_ids_for_name[name]
            else:
                continue
            if mid not in ordered_ids:
                ordered_ids.append(mid)

        if not ordered_ids:
            return pandas.DataFrame(columns=STAR_COLUMNS)

        present_ids = [i for i in ordered_ids if i in set(self.starcache["main_id"])]
        return (
            self.starcache.set_index("main_id", drop=False)
            .reindex(present_ids)
            .reset_index(drop=True)
        )

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
