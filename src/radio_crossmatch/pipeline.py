"""Run a set of catalog matchers over a source list and write the result."""

from __future__ import annotations

import logging
from typing import Dict, Iterable, List, Optional, Sequence

import pandas as pd

from .catalog import CatalogMatcher
from .sources import Source

log = logging.getLogger(__name__)

BASE_COLUMNS = {
    "Name": "Source name from the input list",
    "RA_deg": "Right Ascension (deg, ICRS) used for matching",
    "Dec_deg": "Declination (deg, ICRS) used for matching",
}


def crossmatch(
    sources: Sequence[Source],
    matchers: Iterable[CatalogMatcher],
    on_error: str = "raise",
) -> pd.DataFrame:
    """Cross-match ``sources`` against every matcher.

    Returns a DataFrame with one row per source: ``Name``, any extra input
    columns, ``RA_deg``, ``Dec_deg``, then each matcher's columns in order.
    """
    matchers = list(matchers)
    base = pd.DataFrame(
        [
            {"Name": s.name, **s.extra, "RA_deg": s.ra_deg, "Dec_deg": s.dec_deg}
            for s in sources
        ]
    )
    parts = [base]
    seen = set(base.columns)
    for m in matchers:
        log.info("Matching %s (%d sources)...", m.label, len(sources))
        df = m.match_all(sources, on_error=on_error)
        clash = seen & set(df.columns)
        if clash:
            raise ValueError(
                f"{m.label} would overwrite existing columns {sorted(clash)}; "
                "give it a different prefix/label"
            )
        seen |= set(df.columns)
        parts.append(df.reset_index(drop=True))
    return pd.concat(parts, axis=1)


def column_descriptions(
    matchers: Iterable[CatalogMatcher], extra: Optional[Dict[str, str]] = None
) -> Dict[str, str]:
    """Ordered column descriptions for the output header."""
    desc = dict(BASE_COLUMNS)
    if extra:
        desc.update(extra)
    for m in matchers:
        desc.update(m.columns)
    return desc


def write_output(
    df: pd.DataFrame,
    path: str,
    descriptions: Optional[Dict[str, str]] = None,
    title: str = "Cross-matched source list",
    notes: Sequence[str] = (),
) -> None:
    """Write ``df`` as CSV preceded by a '#'-commented header.

    The header lists each column with its description (where known).  Read
    it back with ``pd.read_csv(path, comment='#')``.
    """
    descriptions = descriptions or {}
    lines = [f"# {title}", "#"]
    lines += [f"# {n}" for n in notes]
    if notes:
        lines.append("#")
    lines += ["# Column descriptions:", "# --------------------"]
    for i, col in enumerate(df.columns, 1):
        lines.append(f"# {i:3d} {col}: {descriptions.get(col, '')}".rstrip())
    lines.append("#")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    df.to_csv(path, mode="a", index=False)
    log.info("Wrote %d rows to %s", len(df), path)
