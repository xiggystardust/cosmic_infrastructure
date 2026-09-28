"""radio_crossmatch: modular cross-matching of source lists against radio
and binary supermassive black hole catalogs."""

from .catalog import (
    CatalogMatcher,
    LocalPositionalCatalog,
    available_catalogs,
    get_catalog,
    register_catalog,
)
from .config import build_matchers, load_config, run_config
from .pipeline import column_descriptions, crossmatch, write_output
from .sources import Source, read_sources, sources_from_dataframe

__version__ = "0.1.0"

__all__ = [
    "CatalogMatcher",
    "LocalPositionalCatalog",
    "Source",
    "available_catalogs",
    "build_matchers",
    "column_descriptions",
    "crossmatch",
    "get_catalog",
    "load_config",
    "read_sources",
    "register_catalog",
    "run_config",
    "sources_from_dataframe",
    "write_output",
]
