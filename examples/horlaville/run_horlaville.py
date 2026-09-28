#!/usr/bin/env python3
"""
Horlaville et al. sample cross-match -- the analysis that radio_crossmatch
grew out of (originally ``full_cross_match.py``, first drafted with Claude
Sonnet 4.5 on 20 July 2026 and then debugged/modified by hand by S. Spolaor).

This shows how to layer a project-specific analysis on top of the generic
package:

1. read a custom input format (the fixed-width MRT table) and attach
   positions from a separate name/RA/Dec file;
2. run the generic cross-match using the catalogs listed in config.yaml;
3. drop known dual/recoiling AGN and split the output by LDA score and by
   presence in the Nyland radio samples.

Usage:  python run_horlaville.py config.yaml
"""

import logging
import sys

import pandas as pd

from radio_crossmatch import (
    build_matchers,
    column_descriptions,
    crossmatch,
    load_config,
    sources_from_dataframe,
    write_output,
)
from radio_crossmatch.config import resolve_path

log = logging.getLogger("horlaville")

# (column, start, end, type) for the fixed-width MRT input
MRT_FIELDS = [
    ("Name", 0, 13, str),
    ("Rank", 14, 17, int),
    ("Distance_Mpc", 18, 23, float),
    ("Survey", 24, 31, str),
    ("logMBH", 32, 37, float),
    ("logh0", 38, 44, float),
    ("LDA", 45, 50, float),
    ("logh0norm", 51, 55, float),
    ("LDAnorm", 56, 60, float),
    ("Score", 61, 65, float),
]

INPUT_DESCRIPTIONS = {
    "Rank": "Total score rank from input catalog",
    "Distance_Mpc": "Luminosity distance in Mpc",
    "Survey": "IFU survey name (MASSIVE or ATLAS3D)",
    "logMBH": "Log of black hole mass in solar masses",
    "logh0": "Log of hypothetical gravitational wave strain",
    "LDA": "Linear Discriminant Analysis score (Equation 2)",
    "logh0norm": "Normalized log of hypothetical GW strain",
    "LDAnorm": "Normalized LDA score",
    "Score": "Total score (Equation 4)",
    "Profile": "Inner light profile (C=core, P=power-law, I=intermediate)",
    "RA_deg": "Right Ascension (deg, J2000) from nameradec_fullsample.txt",
    "Dec_deg": "Declination (deg, J2000) from nameradec_fullsample.txt",
}


def read_mrt(filename):
    """Parse the Horlaville fixed-width table.  Fails loudly on a bad line."""
    rows = []
    with open(filename) as f:
        for n, line in enumerate(f, 1):
            if line.startswith("#") or not line.strip():
                continue
            try:
                row = {c: t(line[a:b].strip()) for c, a, b, t in MRT_FIELDS}
            except ValueError as err:
                raise ValueError(f"{filename}:{n}: could not parse {line!r}: {err}") from err
            row["Profile"] = line[66:67].strip() if len(line) > 66 else ""
            rows.append(row)
    return pd.DataFrame(rows)


def read_radec(filename):
    """name ra dec (deg), whitespace-separated, no header."""
    return pd.read_csv(filename, comment="#", sep=r"\s+", header=None,
                       names=["Name", "ra", "dec"], dtype={"Name": str})


def main(config_path):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    cfg = load_config(config_path)
    base = cfg["_base_dir"]
    hcfg = cfg["horlaville"]

    # 1. inputs ------------------------------------------------------------
    sample = read_mrt(resolve_path(hcfg["sample_file"], base))
    radec = read_radec(resolve_path(hcfg["radec_file"], base))
    missing = set(sample["Name"]) - set(radec["Name"])
    if missing:
        raise SystemExit(f"No RA/Dec in {hcfg['radec_file']} for: {sorted(missing)}")
    merged = sample.merge(radec, on="Name", how="left", validate="one_to_one")
    sources = sources_from_dataframe(merged, name_col="Name", ra_col="ra", dec_col="dec")
    log.info("Read %d galaxies", len(sources))

    # 2. cross-match ---------------------------------------------------------
    matchers = build_matchers(cfg)
    df = crossmatch(sources, matchers, on_error=cfg.get("on_error", "raise"))

    # 3. exclusions and splits ----------------------------------------------
    notes = [f"Match radius: {cfg.get('match_radius_arcsec', 10)} arcsec"]
    for ex in hcfg.get("exclude", []) or []:
        if ex["name"] in set(df["Name"]):
            log.info("Excluding %s: %s", ex["name"], ex.get("reason", ""))
            df = df[df["Name"] != ex["name"]]
            notes.append(f"Excluded {ex['name']}: {ex.get('reason', '')}")
        else:
            log.warning("Exclusion %s not in sample", ex["name"])
    df = df.reset_index(drop=True)

    nyland = (df["NylandMASSIVE"] == "Y") | (df["NylandATLAS"] == "Y")
    lda_gt2 = df["LDA"] > 2
    subsets = {
        "full": (df, "Full sample"),
        "lda_gt2": (df[lda_gt2], "LDA > 2"),
        "lda_lt_m1": (df[df["LDA"] < -1], "LDA < -1"),
        "lda_gt2_ku_only": (df[lda_gt2 & nyland], "LDA > 2, in a Nyland sample (Ku only)"),
        "lda_gt2_c_ku": (df[lda_gt2 & ~nyland], "LDA > 2, not in a Nyland sample (C+Ku)"),
    }

    desc = column_descriptions(matchers, INPUT_DESCRIPTIONS)
    for key, (sub, label) in subsets.items():
        out = resolve_path(hcfg["outputs"][key], base)
        write_output(sub.reset_index(drop=True), out, desc,
                     title=f"Horlaville sample cross-match: {label}", notes=notes)
        print(f"{label:45s} {len(sub):4d} rows -> {out}")

    print(f"LDA > 0: {(df['LDA'] > 0).sum()}   LDA > 2: {lda_gt2.sum()}   "
          f"LDA < -1: {(df['LDA'] < -1).sum()}")
    print(f"In NylandMASSIVE: {(df['NylandMASSIVE'] == 'Y').sum()}   "
          f"In NylandATLAS: {(df['NylandATLAS'] == 'Y').sum()}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config.yaml")
