"""BOBCat: catalog of binary supermassive black hole candidates.

Expects a CSV with (at least) ``name``, ``ra``, ``dec`` columns in degrees
(e.g. ``bobcat_pos_z_uniq.csv``).  Matched by name, then position.
"""

import pandas as pd

from ..catalog import LocalPositionalCatalog


class BOBCat(LocalPositionalCatalog):
    name = "bobcat"
    description = "BOBCat binary SMBH candidate catalog (local CSV: name, ra, dec)"
    columns = {"BOBCat_Match": "Y/N: source is in BOBCat"}
    ra_col, dec_col = "ra", "dec"
    name_cols = ("name",)

    def read_table(self, file):
        return pd.read_csv(file, comment="#")

    def empty_result(self):
        return {"BOBCat_Match": "N"}

    def result_from_row(self, row):
        return {"BOBCat_Match": "Y"}
