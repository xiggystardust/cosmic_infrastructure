"""Build matchers and run a cross-match from a YAML config.

Example::

    match_radius_arcsec: 10
    on_error: raise            # or "warn"
    input:
      file: my_sources.csv     # columns: name, ra, dec (deg)
    output:
      file: crossmatched.csv
    catalogs:
      - type: rfc
        file: data/rfc_2026b_cat.txt
      - type: nvss
        radius_arcsec: 20      # per-catalog override

Relative ``file`` paths are resolved against ``data_dir`` if the config sets
one (itself relative to the config file), otherwise against the config
file's own folder.
"""

from __future__ import annotations

import copy
import logging
import os
from typing import List

import yaml
from astropy import units as u

from .catalog import CatalogMatcher, get_catalog
from .pipeline import column_descriptions, crossmatch, write_output
from .sources import read_sources

log = logging.getLogger(__name__)


def resolve_path(path, base_dir):
    if path is None or os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(base_dir, os.path.expanduser(path)))


def load_config(path: str) -> dict:
    with open(path) as f:
        cfg = yaml.safe_load(f) or {}
    cfg_dir = os.path.dirname(os.path.abspath(path))
    # relative paths resolve against data_dir if given, else the config's folder
    data_dir = cfg.get("data_dir")
    cfg["_base_dir"] = resolve_path(os.path.expanduser(data_dir), cfg_dir) if data_dir else cfg_dir
    return cfg


def build_matchers(cfg: dict) -> List[CatalogMatcher]:
    base_dir = cfg.get("_base_dir", ".")
    default_radius = cfg.get("match_radius_arcsec", 10.0)
    matchers = []
    for entry in cfg.get("catalogs", []):
        entry = copy.deepcopy(entry)
        if entry.pop("enabled", True) is False:
            continue
        cls = get_catalog(entry.pop("type"))
        radius = entry.pop("radius_arcsec", default_radius) * u.arcsec
        if "file" in entry:
            entry["file"] = resolve_path(entry["file"], base_dir)
        matchers.append(cls(radius=radius, **entry))
    return matchers


def run_config(cfg: dict):
    """Run the full cross-match described by ``cfg``; returns the DataFrame."""
    base_dir = cfg.get("_base_dir", ".")
    inp = dict(cfg["input"])
    inp["file"] = resolve_path(inp["file"], base_dir)
    sources = read_sources(**inp)
    log.info("Read %d sources from %s", len(sources), inp["file"])

    matchers = build_matchers(cfg)
    df = crossmatch(sources, matchers, on_error=cfg.get("on_error", "raise"))

    out = cfg.get("output", {})
    if out.get("file"):
        write_output(
            df,
            resolve_path(out["file"], base_dir),
            column_descriptions(matchers, out.get("column_descriptions")),
            title=out.get("title", "Cross-matched source list"),
            notes=[f"Match radius (default): {cfg.get('match_radius_arcsec', 10.0)} arcsec"],
        )
    return df
