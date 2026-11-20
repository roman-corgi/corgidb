import os
import urllib.parse
from typing import Iterable, List, Optional, Union

import numpy as np
import pandas
from astroquery.simbad import SimbadClass

from corgidb.mirrors import SIMBAD_ERRORS, get_simbad
from corgidb.scripts.columns import STAR_COLUMNS
from corgidb.scripts.query_refstars import query_refstars
from corgidb.scripts.query_star import query_star
from corgidb.scripts.query_stars import query_stars
from corgidb.scripts.resolve_star_name import resolve_star_name
from corgidb.utils import get_cache_dir

_DEFAULT_BASE_URL = "https://corgidb.sioslab.com"
_SIMBAD_FIELDS = ("pmra", "pmdec", "plx_value", "rvz_radvel")


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
                deduped = self.starcache.drop_duplicates(
                    subset="main_id", keep="last"
                ).reset_index(drop=True)
                if len(deduped) != len(self.starcache):
                    _atomic_write_csv(deduped, self._starcache_path)
                self.starcache = deduped
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

        fetched = fetched.drop_duplicates(subset="main_id", keep="last")
        self.starcache = pandas.concat(
            [
                self.starcache[~self.starcache["main_id"].isin(fetched["main_id"])],
                fetched,
            ],
            ignore_index=True,
        )
        _atomic_write_csv(self.starcache, self._starcache_path)

    def query_simbad(
        self, st_name: str, simbad: Optional[SimbadClass] = None
    ) -> pandas.DataFrame:
        """Query SIMBAD directly for a star not found in corgidb.

        Used by :meth:`cache_stars` as a fallback for names that don't
        resolve via corgidb's StarAliases table. Only position, proper
        motion, parallax, and radial velocity are available from this
        query; spectype, sy_vmag, sy_imag, and sy_dist are set to null.

        Args:
            st_name (str):
                Star name to look up in SIMBAD.
            simbad (Optional[astroquery.simbad.SimbadClass]):
                Simbad instance to query, as returned by
                :func:`corgidb.mirrors.get_simbad` with the proper motion,
                parallax and radial velocity fields. If None (default), a
                reachable mirror is selected for this call.

        Returns:
            pandas.DataFrame:
                Single-row result with the same columns as query_stars.
                Empty DataFrame if SIMBAD has no match or the query fails.
                If no SIMBAD mirror is reachable, a single row with
                ``st_name`` set and every other column (including
                ``main_id``) NaN.

        """

        unavailable = pandas.DataFrame(
            [{"st_name": st_name}], columns=STAR_COLUMNS, dtype=object
        )

        if simbad is None:
            simbad = get_simbad(*_SIMBAD_FIELDS)
            if simbad is None:
                return unavailable

        try:
            res = simbad.query_object(st_name)
        except SIMBAD_ERRORS:
            print(f"SIMBAD mirror {simbad.server} failed while querying {st_name}.")
            return unavailable
        except Exception:
            res = None

        if res is None or len(res) == 0:
            return pandas.DataFrame(columns=STAR_COLUMNS)

        # Mirrors differ in how they report missing values: CDS returns masked
        # columns, the CfA mirror plain columns holding NaN. Handle both.
        row = res[0]
        missing = {
            col: np.ma.is_masked(row[col]) or np.isnan(row[col])
            for col in ["pmra", "pmdec", "rvz_radvel", "plx_value"]
        }
        plx = None if missing["plx_value"] else float(row["plx_value"])

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
                    "sy_pmra": 0.0 if missing["pmra"] else float(row["pmra"]),
                    "sy_pmdec": 0.0 if missing["pmdec"] else float(row["pmdec"]),
                    "st_radv": 0.0
                    if missing["rvz_radvel"]
                    else float(row["rvz_radvel"]),
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
        collapse to a single row. If no SIMBAD mirror is reachable, each
        name that needed SIMBAD gets a row with ``st_name`` set and every
        other column NaN; these rows are never written to the cache.

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
        # all-NaN rows for names SIMBAD could not be asked about; never cached
        unavailable_rows: dict[str, pandas.DataFrame] = {}
        newly_fetched_ids: set = set()
        new_simbad_rows = []
        simbad = None
        simbad_checked = False
        for name in st_names:
            if name in resolved or name in unavailable_rows:
                continue
            if not force_new and name in cached_ids:
                # already cached from an earlier SIMBAD lookup under this
                # exact name
                simbad_ids_for_name[name] = name
                continue
            if not simbad_checked:
                # pick a mirror once per call, and only if SIMBAD is needed
                simbad = get_simbad(*_SIMBAD_FIELDS)
                simbad_checked = True
            if simbad is None:
                unavailable_rows[name] = pandas.DataFrame(
                    [{"st_name": name}], columns=STAR_COLUMNS, dtype=object
                )
                continue
            simbad_row = self.query_simbad(name, simbad=simbad)
            if simbad_row.empty:
                continue
            mid = simbad_row.iloc[0]["main_id"]
            if pandas.isna(mid):
                unavailable_rows[name] = simbad_row
                continue
            simbad_ids_for_name[name] = mid
            # a different name earlier in this same batch may have already
            # resolved to the same SIMBAD main_id; only queue one row for it
            if mid not in cached_ids and mid not in newly_fetched_ids:
                newly_fetched_ids.add(mid)
                new_simbad_rows.append(simbad_row)

        if new_simbad_rows:
            self._merge_into_cache(pandas.concat(new_simbad_rows, ignore_index=True))

        starcache_unique = self.starcache.drop_duplicates(
            subset="main_id", keep="last"
        ).set_index("main_id", drop=False)
        seen_ids: set = set()
        pieces: List[pandas.DataFrame] = []
        for name in st_names:
            if name in unavailable_rows:
                pieces.append(unavailable_rows.pop(name))
                continue
            if name in resolved:
                mid = resolved[name]
            elif name in simbad_ids_for_name:
                mid = simbad_ids_for_name[name]
            else:
                continue
            if mid in seen_ids or mid not in starcache_unique.index:
                continue
            seen_ids.add(mid)
            pieces.append(starcache_unique.loc[[mid]])

        if not pieces:
            return pandas.DataFrame(columns=STAR_COLUMNS)

        return pandas.concat(pieces, ignore_index=True)[STAR_COLUMNS]

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
