"""TEMPLATE for a new catalog module -- copy to ``mycatalog.py`` and edit.

Files starting with an underscore are not auto-imported, so this template is
never registered.  Once copied (without the underscore) and given a ``name``,
your class is available as ``type: mycatalog`` in a config file.

Pick the base class that fits:

* LocalPositionalCatalog -- a catalog file on disk with RA/Dec (and maybe
  names).  Implement read_table() and result_from_row().  Example below.
* VizierCatalog (catalogs/vizier.py) -- a VizieR table.  Usually just class
  attributes; or use ``type: vizier`` in the config with no code.
* CatalogMatcher -- anything else (a TAP service, a web API, image cutouts).
  Implement match(source) -> dict.
"""

import numpy as np
import pandas as pd
from astropy import units as u

from ..catalog import LocalPositionalCatalog


class MyCatalog(LocalPositionalCatalog):
    # Registry key used in configs.  Must be unique.
    name = "mycatalog"
    description = "One line describing the catalog and where to get it"

    # Output columns, in order, with descriptions for the output header.
    # Prefix them with the catalog name so they cannot collide.
    columns = {
        "MYCAT_Match": "Y/N: source is in MyCatalog",
        "MYCAT_Flux_mJy": "MyCatalog 5 GHz flux density (mJy)",
    }

    # Column names/units of the positions in the table read_table() returns.
    ra_col, dec_col = "RAJ2000", "DEJ2000"
    coord_unit = (u.deg, u.deg)       # (u.hourangle, u.deg) for hh:mm:ss
    name_cols = ()                    # e.g. ("Name",) to also match by name

    def read_table(self, file):
        return pd.read_csv(file, comment="#")

    def empty_result(self):
        return {"MYCAT_Match": "N", "MYCAT_Flux_mJy": np.nan}

    def result_from_row(self, row):
        return {"MYCAT_Match": "Y", "MYCAT_Flux_mJy": row["S5GHz"]}
