"""RFC: Radio Fundamental Catalog (VLBI), http://astrogeo.org/rfc/

Expects the fixed-width ``rfc_XXXXx_cat.txt`` file.  Position-only match;
returns short- (Fs) and long-baseline (Fl) flux densities in Jy for S, C, X,
U and K bands.  Values <= -9 in the file (no data) become NaN.
"""

import logging

import numpy as np
import pandas as pd

from ..catalog import LocalPositionalCatalog

log = logging.getLogger(__name__)

# (column suffix, start, end) of each flux field in the fixed-width file
_FLUX_FIELDS = [
    ("FsS", 104, 110), ("FlS", 118, 124),
    ("FsC", 126, 132), ("FlC", 140, 146),
    ("FsX", 148, 154), ("FlX", 162, 168),
    ("FsU", 170, 176), ("FlU", 184, 190),
    ("FsK", 192, 198), ("FlK", 206, 212),
]


def _parse_flux(s):
    try:
        val = float(s)
    except ValueError:
        return np.nan
    return val if val > -9.0 else np.nan


def _flux_desc(key):
    base = "short" if key[1] == "s" else "long"
    return f"RFC {key[2]}-band {base}-baseline flux density (Jy)"


class RFC(LocalPositionalCatalog):
    name = "rfc"
    description = "Radio Fundamental Catalog, VLBI flux densities (local fixed-width file)"
    columns = {"RFC_Match": "Y/N: source is in the Radio Fundamental Catalog"}
    columns.update({f"RFC_{k}": _flux_desc(k) for k, _, _ in _FLUX_FIELDS})
    ra_col, dec_col = "RFC_RA", "RFC_Dec"

    def __init__(self, file, **kwargs):
        kwargs.setdefault("name_match", "none")
        super().__init__(file, **kwargs)

    def read_table(self, file):
        rows, bad = [], 0
        with open(file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                try:
                    ra_h = int(line[26:28]) + int(line[29:31]) / 60 + float(line[32:41]) / 3600
                    sign = -1.0 if line[42:43] == "-" else 1.0
                    dec_d = int(line[43:45]) + int(line[46:48]) / 60 + float(line[49:57]) / 3600
                    row = {
                        "RFC_Name": line[0:14].strip(),
                        "RFC_CommonName": line[16:24].strip(),
                        "RFC_RA": 15.0 * ra_h,
                        "RFC_Dec": sign * dec_d,
                    }
                    for key, a, b in _FLUX_FIELDS:
                        row[f"RFC_{key}"] = _parse_flux(line[a:b]) if len(line) > b else np.nan
                    rows.append(row)
                except Exception:  # noqa: BLE001
                    bad += 1
        if bad:
            log.warning("RFC: skipped %d unparseable lines in %s", bad, file)
        return pd.DataFrame(rows)

    def empty_result(self):
        res = {c: np.nan for c in self.columns}
        res["RFC_Match"] = "N"
        return res

    def result_from_row(self, row):
        res = {c: row[c] for c in self.columns if c != "RFC_Match"}
        res["RFC_Match"] = "Y"
        return {c: res[c] for c in self.columns}
