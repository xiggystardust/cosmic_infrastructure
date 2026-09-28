"""Small shared helpers."""

import re

_WS = re.compile(r"\s+")


def normalize_name(name) -> str:
    """Upper-case and strip all whitespace, so 'NGC 4486' == 'ngc4486'."""
    return _WS.sub("", str(name)).upper()


def yn(flag: bool) -> str:
    """Convert a boolean to the 'Y'/'N' convention used in output tables."""
    return "Y" if flag else "N"
