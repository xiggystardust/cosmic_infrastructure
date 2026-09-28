"""Catalogs queried on the fly from VizieR (VLASS, FIRST, NVSS, or any other).

:class:`VizierCatalog` handles the query, picks the *nearest* row within the
match radius (rows are sorted by distance server-side via ``+_r``), and
optionally reports whether the position lies in the survey footprint.

Adding a new VizieR survey usually takes no code at all -- use the generic
``vizier`` type in the config::

    - type: vizier
      catalog: VIII/81B/sumss212     # VizieR table ID
      prefix: SUMSS
      flux_column: St              # column to report as <prefix>_Flux_mJy
      dec_max: -30                   # footprint (optional)
      extra_columns: [MajAxis]       # copied as <prefix>_<col>

or subclass :class:`VizierCatalog` with class attributes, as below.
"""

from __future__ import annotations

import logging
import time
from typing import Optional, Sequence

import numpy as np

from ..catalog import CatalogMatcher
from ..sources import Source
from ..utils import yn

log = logging.getLogger(__name__)


class VizierCatalog(CatalogMatcher):
    """Base class for a single VizieR table.

    Class attributes (all overridable as constructor keyword arguments):

    vizier_id      VizieR table identifier, e.g. "VIII/65/nvss"
    prefix         Output column prefix, e.g. "NVSS"
    flux_column    VizieR column reported as ``<prefix>_Flux_mJy`` (or None)
    flux_desc      Description of that column for the output header
    dec_min/max    Approximate declination limits of the footprint, or None
    report_footprint  Add a ``<prefix>_InFootprint`` Y/N column
    nondetection_value  Flux reported for an in-footprint non-detection
    extra_columns  More VizieR columns to copy as ``<prefix>_<col>``
    """

    vizier_id: str = ""
    prefix: str = ""
    flux_column: Optional[str] = None
    flux_desc: str = "flux density (mJy)"
    dec_min: Optional[float] = None
    dec_max: Optional[float] = None
    report_footprint: bool = True
    nondetection_value: float = -999.0
    extra_columns: Sequence[str] = ()
    retries: int = 2

    _OPTIONS = (
        "vizier_id", "prefix", "flux_column", "flux_desc", "dec_min", "dec_max",
        "report_footprint", "nondetection_value", "extra_columns", "retries",
    )

    def __init__(self, **kwargs):
        # accept "catalog" as an alias for vizier_id in configs
        if "catalog" in kwargs:
            kwargs["vizier_id"] = kwargs.pop("catalog")
        for opt in self._OPTIONS:
            if opt in kwargs:
                setattr(self, opt, kwargs.pop(opt))
        super().__init__(**kwargs)
        if not self.vizier_id or not self.prefix:
            raise ValueError(f"{type(self).__name__} needs vizier_id and prefix")
        self.columns = self._build_columns()
        self._vizier = None

    # columns depend on options, so build per instance
    def _build_columns(self):
        p = self.prefix
        cols = {}
        if self.report_footprint:
            lims = []
            if self.dec_min is not None:
                lims.append(f"dec > {self.dec_min:g}")
            if self.dec_max is not None:
                lims.append(f"dec < {self.dec_max:g}")
            where = f" ({', '.join(lims)})" if lims else ""
            cols[f"{p}_InFootprint"] = (
                f"Y/N: inside the approximate {p} footprint{where} or detected"
            )
        if self.flux_column:
            nd = (
                f"; {self.nondetection_value:g} = non-detection in footprint"
                if self.report_footprint else ""
            )
            cols[f"{p}_Flux_mJy"] = f"{p} {self.flux_desc} [{self.flux_column}]{nd}"
        for c in self.extra_columns:
            cols[f"{p}_{c}"] = f"{p} catalog column {c}"
        return cols

    @property
    def label(self):
        return self.prefix

    def load(self):
        from astroquery.vizier import Vizier

        self._vizier = Vizier(columns=["**", "+_r"], row_limit=1)
        # A wrong table ID makes query_region return nothing, which would
        # look exactly like "no detection" for every source.  Check it once.
        found = Vizier(row_limit=1).get_catalogs(self.vizier_id)
        if self.vizier_id not in found.keys():
            raise ValueError(
                f"VizieR table '{self.vizier_id}' not found for {self.prefix} "
                f"(VizieR returned {list(found.keys()) or 'nothing'}). "
                "Check the full table ID, e.g. 'VIII/92/first14'."
            )

    def in_footprint(self, coord) -> bool:
        dec = coord.dec.deg
        if self.dec_min is not None and dec < self.dec_min:
            return False
        if self.dec_max is not None and dec > self.dec_max:
            return False
        return True

    def query(self, coord):
        """Return the nearest VizieR row within the radius, or None."""
        for attempt in range(self.retries + 1):
            try:
                res = self._vizier.query_region(coord, radius=self.radius, catalog=self.vizier_id)
                break
            except Exception as err:  # noqa: BLE001
                if attempt == self.retries:
                    raise
                log.warning("%s query failed (%s); retrying", self.prefix, err)
                time.sleep(2 * (attempt + 1))
        if len(res) == 0 or len(res[0]) == 0:
            return None
        return res[0][0]

    def match(self, source: Source) -> dict:
        self.ensure_loaded()
        p = self.prefix
        result = self.empty_result()
        if source.coord is None:
            if self.report_footprint:
                result[f"{p}_InFootprint"] = "N"
            return result

        row = self.query(source.coord)
        inside = self.in_footprint(source.coord)
        if self.report_footprint:
            # a detection implies coverage, whatever the approximate limits say
            result[f"{p}_InFootprint"] = yn(inside or row is not None)

        if row is not None:
            if self.flux_column:
                colnames = row.colnames
                result[f"{p}_Flux_mJy"] = (
                    float(row[self.flux_column]) if self.flux_column in colnames else np.nan
                )
            for c in self.extra_columns:
                result[f"{p}_{c}"] = row[c] if c in row.colnames else np.nan
        elif inside and self.report_footprint and self.flux_column:
            result[f"{p}_Flux_mJy"] = self.nondetection_value
        return result


class GenericVizier(VizierCatalog):
    """Any VizieR table, configured entirely from the config file."""

    name = "vizier"
    description = "Any VizieR table (set catalog, prefix, flux_column, dec_min/max in config)"


class VLASS(VizierCatalog):
    name = "vlass"
    description = "VLASS epoch 1 Quick Look components, Gordon+2021 (VizieR J/ApJS/255/30/comp)"
    vizier_id = "J/ApJS/255/30/comp"
    prefix = "VLASS"
    flux_column = "Ftot"
    flux_desc = "3 GHz total flux density (mJy)"
    dec_min = -40.0
    report_footprint = False  # matches the original script's output


class FIRST(VizierCatalog):
    name = "first"
    description = "FIRST 1.4 GHz survey, 14Dec17 catalog (VizieR VIII/92/first14)"
    vizier_id = "VIII/92/first14"
    prefix = "FIRST"
    flux_column = "Fpeak"
    flux_desc = "1.4 GHz peak flux density (mJy/beam)"
    # Crude: FIRST covers ~10,575 deg^2, mostly the north Galactic cap and an
    # equatorial strip, not the whole sky north of -10 deg.
    dec_min = -10.0


class NVSS(VizierCatalog):
    name = "nvss"
    description = "NVSS 1.4 GHz survey (VizieR VIII/65)"
    vizier_id = "VIII/65/nvss"
    prefix = "NVSS"
    flux_column = "S1.4"
    flux_desc = "1.4 GHz integrated flux density (mJy)"
    dec_min = -40.0
