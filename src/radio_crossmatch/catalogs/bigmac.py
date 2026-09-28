"""BIG MAC: Big Multi-AGN Catalog of binary/dual AGN candidates.

Expects the BIG MAC main-table CSV (e.g. ``BigMAC_maintable_DR0p9.csv``),
whose first line is a title row (hence ``header=1``).  Sources are matched by
name against Name1/Name2, then by position against RA1/Dec1 (sexagesimal).
"""

import pandas as pd
from astropy import units as u

from ..catalog import LocalPositionalCatalog


class BigMAC(LocalPositionalCatalog):
    name = "bigmac"
    description = "BIG MAC binary/dual AGN candidate catalog (local CSV)"
    columns = {
        "BIGMAC_detected": "Y/N: source is in BIG MAC",
        "BIGMAC_PrimaryType": "Primary system type from BIG MAC",
        "BIGMAC_SecondaryType": "Secondary system type from BIG MAC",
        "BIGMAC_TertiaryType": "Tertiary system type from BIG MAC",
    }
    ra_col, dec_col = "RA1", "Dec1"
    coord_unit = (u.hourangle, u.deg)
    name_cols = ("Name1", "Name2")

    def __init__(self, file, header_row: int = 1, **kwargs):
        super().__init__(file, **kwargs)
        self.header_row = header_row

    def read_table(self, file):
        return pd.read_csv(file, header=self.header_row)

    def empty_result(self):
        return {
            "BIGMAC_detected": "N",
            "BIGMAC_PrimaryType": "",
            "BIGMAC_SecondaryType": "",
            "BIGMAC_TertiaryType": "",
        }

    def result_from_row(self, row):
        return {
            "BIGMAC_detected": "Y",
            "BIGMAC_PrimaryType": row.get("Primary System Type", ""),
            "BIGMAC_SecondaryType": row.get("Secondary System Type", ""),
            "BIGMAC_TertiaryType": row.get("Tertiary System Type", ""),
        }
