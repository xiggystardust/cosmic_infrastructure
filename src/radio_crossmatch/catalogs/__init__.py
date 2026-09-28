"""Built-in catalog modules.

Every module in this package whose name does not start with an underscore is
imported automatically, which registers the CatalogMatcher subclasses it
defines.  To add a catalog, drop a new file here (see ``_template.py``).
"""

import importlib
import pkgutil

for _mod in pkgutil.iter_modules(__path__):
    if not _mod.name.startswith("_"):
        importlib.import_module(f"{__name__}.{_mod.name}")
