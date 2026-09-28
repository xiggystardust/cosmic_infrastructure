"""Command-line entry point: ``radio-crossmatch config.yaml``."""

import argparse
import logging
import sys

from .catalog import available_catalogs, get_catalog
from .config import load_config, run_config


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="radio-crossmatch",
        description="Cross-match a source list against radio and binary-SMBH catalogs.",
    )
    p.add_argument("config", nargs="?", help="YAML config file")
    p.add_argument("--list", action="store_true", help="list available catalog modules")
    p.add_argument("-v", "--verbose", action="count", default=0)
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose > 1 else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    if args.list:
        for name in available_catalogs():
            cls = get_catalog(name)
            print(f"{name:12s} {cls.description}")
        return 0
    if not args.config:
        p.print_help()
        return 1
    run_config(load_config(args.config))
    return 0


if __name__ == "__main__":
    sys.exit(main())
