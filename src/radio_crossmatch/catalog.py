"""
Base classes and registry for catalog cross-match modules.

Every catalog is a subclass of :class:`CatalogMatcher`.  Giving the subclass a
``name`` class attribute registers it automatically, so it can be selected by
that name in a YAML config (``type: <name>``) or via :func:`get_catalog`.

Two helper base classes cover the common cases:

* :class:`LocalPositionalCatalog` -- a catalog file on disk with RA/Dec
  columns (optionally names).  Subclasses only implement ``read_table``.
* :class:`~radio_crossmatch.catalogs.vizier.VizierCatalog` -- a catalog
  queried on the fly from VizieR.

See CONTRIBUTING.md for a walkthrough of writing a new module.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Sequence, Type

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from .sources import Source
from .utils import normalize_name

log = logging.getLogger(__name__)

_REGISTRY: Dict[str, Type["CatalogMatcher"]] = {}


def register_catalog(cls: Type["CatalogMatcher"]) -> Type["CatalogMatcher"]:
    """Register a CatalogMatcher subclass under ``cls.name``.

    Called automatically for any subclass that defines ``name``; can also be
    used as a decorator for classes defined outside this package.
    """
    key = cls.name.lower()
    if key in _REGISTRY and _REGISTRY[key] is not cls:
        raise ValueError(
            f"Catalog name '{cls.name}' is already registered by "
            f"{_REGISTRY[key].__module__}.{_REGISTRY[key].__name__}"
        )
    _REGISTRY[key] = cls
    return cls


def get_catalog(name: str) -> Type["CatalogMatcher"]:
    """Look up a registered catalog class by name (case-insensitive)."""
    from . import catalogs  # noqa: F401  (ensures built-ins are registered)

    try:
        return _REGISTRY[name.lower()]
    except KeyError:
        raise KeyError(
            f"Unknown catalog '{name}'. Available: {', '.join(available_catalogs())}"
        ) from None


def available_catalogs() -> List[str]:
    """Names of all registered catalogs."""
    from . import catalogs  # noqa: F401

    return sorted(_REGISTRY)


class CatalogMatcher(ABC):
    """Base class for one cross-match module.

    Subclasses set:

    ``name``
        Registry key used in configs (e.g. ``"rfc"``).  Leave unset on
        abstract intermediate classes.
    ``description``
        One-line description of the catalog.
    ``columns``
        Ordered mapping ``{output_column: description}``.  Every dict returned
        by :meth:`match` must contain exactly these keys; the descriptions are
        written into the output file header.

    and implement :meth:`match` (per source) and optionally :meth:`load`
    (one-off setup such as reading a file) or :meth:`match_all` (vectorized
    matching of the whole sample).
    """

    name: Optional[str] = None
    description: str = ""
    columns: Dict[str, str] = {}

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.__dict__.get("name"):
            register_catalog(cls)

    def __init__(self, radius: u.Quantity = 10 * u.arcsec, **options):
        self.radius = u.Quantity(radius, u.arcsec)
        if options:
            raise TypeError(
                f"{type(self).__name__} got unexpected options: {sorted(options)}"
            )
        self._loaded = False

    @property
    def label(self) -> str:
        """Human-readable name used in log messages."""
        return self.name or type(self).__name__

    # -- lifecycle ----------------------------------------------------------

    def load(self) -> None:
        """Do any one-off setup (read files, open connections)."""

    def ensure_loaded(self) -> None:
        if not self._loaded:
            self.load()
            self._loaded = True

    # -- matching -----------------------------------------------------------

    def empty_result(self) -> dict:
        """Result for a non-match.  Override if NaN is not a sensible default."""
        return {col: np.nan for col in self.columns}

    @abstractmethod
    def match(self, source: Source) -> dict:
        """Return ``{column: value}`` for one source (keys == ``self.columns``)."""

    def match_all(self, sources: Sequence[Source], on_error: str = "raise") -> pd.DataFrame:
        """Match every source.  Default loops over :meth:`match`.

        ``on_error="raise"`` stops on the first failure (the original
        script's behaviour); ``"warn"`` logs it and fills that row with
        :meth:`empty_result`.
        """
        self.ensure_loaded()
        rows = []
        for n, src in enumerate(sources, 1):
            log.debug("  %s [%d/%d] %s", self.label, n, len(sources), src.name)
            try:
                res = self.match(src)
            except Exception as err:  # noqa: BLE001
                if on_error == "raise":
                    raise RuntimeError(f"{self.label} failed on {src.name}: {err}") from err
                log.warning("%s failed on %s: %s", self.label, src.name, err)
                res = self.empty_result()
            missing = set(self.columns) - set(res)
            if missing:
                raise KeyError(f"{self.label}.match() did not return columns {sorted(missing)}")
            rows.append(res)
        return pd.DataFrame(rows, columns=list(self.columns))

    def __repr__(self):
        return f"<{type(self).__name__} name={self.name!r} radius={self.radius}>"


class LocalPositionalCatalog(CatalogMatcher):
    """Helper base for catalogs stored in a local file.

    Subclasses implement :meth:`read_table`, returning a DataFrame, and set
    ``ra_col``/``dec_col``/``coord_unit`` (plus optionally ``name_cols``).
    Matching is done by name first (if ``name_cols`` is set and
    ``name_match`` != ``"none"``), then by nearest position within ``radius``.
    Subclasses implement :meth:`result_from_row` to turn the matched catalog
    row into output columns.

    Positions are converted to a SkyCoord once at load time and matched with
    ``match_to_catalog_sky``, so large catalogs are fast.
    """

    ra_col: str = "ra"
    dec_col: str = "dec"
    coord_unit = (u.deg, u.deg)
    name_cols: Sequence[str] = ()

    def __init__(self, file: str, name_match: str = "contains", **kwargs):
        """
        Parameters
        ----------
        file : str
            Path to the catalog file.
        name_match : {"contains", "exact", "none"}
            How input names are compared to catalog names (after removing
            whitespace and upper-casing).  ``"contains"`` reproduces the
            original script: a catalog entry matches if it *contains* the
            input name -- note this means ``NGC410`` matches ``NGC4101``.
        """
        super().__init__(**kwargs)
        if name_match not in ("contains", "exact", "none"):
            raise ValueError("name_match must be 'contains', 'exact' or 'none'")
        self.file = file
        self.name_match = name_match
        self.table: pd.DataFrame = pd.DataFrame()
        self.coords: Optional[SkyCoord] = None
        self._norm_names: Dict[str, pd.Series] = {}

    @abstractmethod
    def read_table(self, file: str) -> pd.DataFrame:
        """Read the catalog file into a DataFrame."""

    @abstractmethod
    def result_from_row(self, row: pd.Series) -> dict:
        """Convert a matched catalog row to the output columns."""

    def load(self):
        self.table = self.read_table(self.file).reset_index(drop=True)
        log.info("Loaded %d rows from %s (%s)", len(self.table), self.file, self.name)
        if self.ra_col in self.table and self.dec_col in self.table and len(self.table):
            self.coords = SkyCoord(
                self.table[self.ra_col].values,
                self.table[self.dec_col].values,
                unit=self.coord_unit,
                frame="icrs",
            )
        for col in self.name_cols:
            if col in self.table:
                self._norm_names[col] = self.table[col].astype(str).map(normalize_name)

    def find_by_name(self, name: str) -> Optional[int]:
        if self.name_match == "none" or not name:
            return None
        target = normalize_name(name)
        for col, names in self._norm_names.items():
            if self.name_match == "exact":
                mask = names == target
            else:
                mask = names.str.contains(target, regex=False)
            if mask.any():
                return int(np.flatnonzero(mask.values)[0])
        return None

    def find_by_position(self, coord: Optional[SkyCoord]) -> Optional[int]:
        if coord is None or self.coords is None:
            return None
        sep = coord.separation(self.coords)
        i = int(sep.argmin())
        return i if sep[i] < self.radius else None

    def match(self, source: Source) -> dict:
        self.ensure_loaded()
        i = self.find_by_name(source.name)
        how = "name"
        if i is None:
            i = self.find_by_position(source.coord)
            how = "position"
        if i is None:
            return self.empty_result()
        log.info("  %s: %s matched by %s", self.label, source.name, how)
        return self.result_from_row(self.table.iloc[i])
