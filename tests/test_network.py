"""Live VizieR checks.  Skip with:  pytest -m 'not network'"""

import pytest
from astropy.coordinates import SkyCoord

import radio_crossmatch as rx

pytestmark = pytest.mark.network

M87 = rx.Source("M87", SkyCoord(187.70593, 12.39112, unit="deg"))


@pytest.mark.parametrize("name", ["nvss", "first", "vlass"])
def test_builtin_vizier_tables_exist_and_detect_m87(name):
    m = rx.get_catalog(name)()
    r = m.match(M87)          # load() validates the table ID
    assert r[f"{m.prefix}_Flux_mJy"] > 0


def test_bad_table_id_is_caught():
    m = rx.get_catalog("vizier")(catalog="VIII/92/first", prefix="BAD", flux_column="Fpeak")
    with pytest.raises(ValueError, match="not found"):
        m.ensure_loaded()
