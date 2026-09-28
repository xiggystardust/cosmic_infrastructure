"""
Input source lists.

The pipeline works on a list of :class:`Source` objects: a name, a sky
position, and any extra per-source columns that should be carried through to
the output unchanged.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

log = logging.getLogger(__name__)


@dataclass
class Source:
    name: str
    coord: Optional[SkyCoord]
    extra: Dict = field(default_factory=dict)

    @property
    def ra_deg(self):
        return self.coord.ra.deg if self.coord is not None else float("nan")

    @property
    def dec_deg(self):
        return self.coord.dec.deg if self.coord is not None else float("nan")


def resolve_name(name: str) -> Optional[SkyCoord]:
    """Resolve a name with CDS Sesame.

    Caution: name resolvers can confuse close pairs (e.g. NGC7436 vs
    NGC7436B), so supplying explicit positions is preferred.
    """
    try:
        return SkyCoord.from_name(name)
    except Exception as err:  # noqa: BLE001
        log.warning("Could not resolve %s: %s", name, err)
        return None


def sources_from_dataframe(
    df: pd.DataFrame,
    name_col: str = "name",
    ra_col: Optional[str] = "ra",
    dec_col: Optional[str] = "dec",
    coord_unit: str = "deg",
    resolve_missing: bool = False,
) -> List[Source]:
    """Build Sources from a DataFrame.

    Parameters
    ----------
    coord_unit : "deg" or "hourangle"
        Unit of the RA column ("hourangle" for sexagesimal hh:mm:ss strings;
        Dec is always degrees / dd:mm:ss).
    resolve_missing : bool
        If True, rows with missing RA/Dec are resolved by name via Sesame.
    """
    ra_unit = {"deg": u.deg, "hourangle": u.hourangle}[coord_unit]
    has_pos = ra_col in df.columns and dec_col in df.columns
    sources = []
    for _, row in df.iterrows():
        name = str(row[name_col]).strip()
        coord = None
        if has_pos and pd.notna(row[ra_col]) and pd.notna(row[dec_col]):
            coord = SkyCoord(row[ra_col], row[dec_col], unit=(ra_unit, u.deg), frame="icrs")
        elif resolve_missing:
            coord = resolve_name(name)
        if coord is None:
            log.warning("No position for %s", name)
        extra = {k: v for k, v in row.items() if k not in (name_col, ra_col, dec_col)}
        sources.append(Source(name, coord, extra))
    return sources


def read_sources(
    file: str,
    name_col: str = "name",
    ra_col: str = "ra",
    dec_col: str = "dec",
    coord_unit: str = "deg",
    format: str = "auto",
    resolve_missing: bool = False,
) -> List[Source]:
    """Read a source list from disk.

    ``format`` is "csv", "whitespace", or "auto" (csv if the file name ends
    in .csv, otherwise whitespace-delimited).  Lines starting with '#' are
    ignored.  For whitespace files with no header row, pass
    ``name_col=0, ra_col=1, dec_col=2`` (integer column positions).
    """
    if format == "auto":
        format = "csv" if str(file).lower().endswith(".csv") else "whitespace"
    header = None if isinstance(name_col, int) else "infer"
    if format == "csv":
        df = pd.read_csv(file, comment="#", header=header)
    elif format == "whitespace":
        df = pd.read_csv(file, comment="#", sep=r"\s+", header=header)
    else:
        raise ValueError(f"Unknown format {format!r}")
    if header is None:  # positional columns -> give them names
        df = df.rename(columns={name_col: "name", ra_col: "ra", dec_col: "dec"})
        df.columns = [str(c) for c in df.columns]
        name_col, ra_col, dec_col = "name", "ra", "dec"
    return sources_from_dataframe(df, name_col, ra_col, dec_col, coord_unit, resolve_missing)
