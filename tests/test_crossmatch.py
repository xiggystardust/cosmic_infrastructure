"""Offline tests: registry, local catalogs, pipeline, config.  No network."""

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

import radio_crossmatch as rx
from radio_crossmatch.catalog import CatalogMatcher
from radio_crossmatch.sources import Source


def src(name, ra, dec):
    return Source(name, SkyCoord(ra * u.deg, dec * u.deg))


def test_builtin_catalogs_registered():
    names = rx.available_catalogs()
    for n in ["bigmac", "bobcat", "rfc", "namelist", "vlass", "first", "nvss", "vizier"]:
        assert n in names
    assert "mycatalog" not in names  # template is not auto-imported


def test_subclass_autoregisters_and_duplicate_rejected():
    class Dummy(CatalogMatcher):
        name = "dummy_test"
        columns = {"D": "dummy"}

        def match(self, source):
            return {"D": source.name}

    assert rx.get_catalog("DUMMY_TEST") is Dummy
    with pytest.raises(ValueError):
        class Dummy2(CatalogMatcher):  # noqa: F841
            name = "dummy_test"

            def match(self, source):
                return {}


def test_bobcat_name_and_position(tmp_path):
    f = tmp_path / "bobcat.csv"
    f.write_text("name,ra,dec\nNGC 1234,10.0,20.0\nFOO,50.0,-5.0\n")
    m = rx.get_catalog("bobcat")(file=str(f))
    df = m.match_all([
        src("ngc1234", 0, 0),                    # name match (whitespace/case)
        src("X", 50.0 + 5 / 3600, -5.0),         # 5" away -> position match
        src("Y", 50.0 + 20 / 3600, -5.0),        # 20" away -> no match
    ])
    assert list(df["BOBCat_Match"]) == ["Y", "Y", "N"]


def test_name_match_modes(tmp_path):
    f = tmp_path / "bobcat.csv"
    f.write_text("name,ra,dec\nNGC4101,10.0,20.0\n")
    s = src("NGC410", 100, 0)
    assert rx.get_catalog("bobcat")(file=str(f)).match(s)["BOBCat_Match"] == "Y"
    exact = rx.get_catalog("bobcat")(file=str(f), name_match="exact")
    assert exact.match(s)["BOBCat_Match"] == "N"


def test_bigmac(tmp_path):
    f = tmp_path / "bigmac.csv"
    f.write_text(
        "BIG MAC DR0.9\n"
        "Name1,Name2,RA1,Dec1,Primary System Type,Secondary System Type,Tertiary System Type\n"
        "OJ 287,,08:54:48.87,+20:06:30.6,Binary,Blazar,\n"
    )
    m = rx.get_catalog("bigmac")(file=str(f))
    c = SkyCoord("08:54:48.87 +20:06:30.6", unit=(u.hourangle, u.deg))
    r = m.match(Source("J0854+2006", c))
    assert r["BIGMAC_detected"] == "Y" and r["BIGMAC_PrimaryType"] == "Binary"
    assert m.match(src("nothing", 0, 0))["BIGMAC_detected"] == "N"


def _rfc_line(name, ra, dec, fs_s="  0.50", fl_c="  0.25"):
    buf = [" "] * 215
    def put(s, a):
        buf[a:a + len(s)] = list(s)
    put(name, 0); put("0851+202", 16)
    put(ra[0], 26); put(ra[1], 29); put(ra[2], 32)
    put(dec[0], 42); put(dec[1], 43); put(dec[2], 46); put(dec[3], 49)
    put(fs_s, 104); put(fl_c, 140); put(" -1.00", 148)
    return "".join(buf) + "\n"


def test_rfc_parse_and_match(tmp_path):
    f = tmp_path / "rfc.txt"
    f.write_text("# header\n" + _rfc_line("J0854+2006", ("08", "54", "48.874930"),
                                        ("+", "20", "06", "30.64089")))
    m = rx.get_catalog("rfc")(file=str(f))
    m.ensure_loaded()
    assert m.table["RFC_RA"][0] == pytest.approx(133.70364554, abs=1e-6)
    assert m.table["RFC_Dec"][0] == pytest.approx(20.10851136, abs=1e-6)
    r = m.match(src("OJ287", 133.70364554, 20.10851136))
    assert r["RFC_Match"] == "Y"
    assert r["RFC_FsS"] == 0.5 and r["RFC_FlC"] == 0.25
    assert r["RFC_FsX"] == -1.0        # only values <= -9 are treated as missing
    assert np.isnan(r["RFC_FlX"])      # blank field
    assert list(r) == list(m.columns)


def test_namelist(tmp_path):
    ws = tmp_path / "atlas.txt"
    ws.write_text("# name other\nNGC4486 1 2\nNGC 4278 3 4\n")
    fx = tmp_path / "fixed.txt"
    fx.write_text("NGC0383       1 2\n")
    a = rx.get_catalog("namelist")(file=str(ws), label="A")
    b = rx.get_catalog("namelist")(file=str(fx), label="B", format="fixed", start=0, end=13)
    assert a.match(src("NGC4486", 0, 0)) == {"A": "Y"}
    assert a.match(src("NGC0383", 0, 0)) == {"A": "N"}
    assert b.match(src("NGC0383", 0, 0)) == {"B": "Y"}


def test_vizier_columns_and_footprint():
    first = rx.get_catalog("first")()
    assert list(first.columns) == ["FIRST_InFootprint", "FIRST_Flux_mJy"]
    assert list(rx.get_catalog("vlass")().columns) == ["VLASS_Flux_mJy"]
    g = rx.get_catalog("vizier")(catalog="J/x/y", prefix="SUMSS", flux_column="St",
                                  dec_max=-30, extra_columns=["MajAxis"])
    assert list(g.columns) == ["SUMSS_InFootprint", "SUMSS_Flux_mJy", "SUMSS_MajAxis"]
    assert g.in_footprint(SkyCoord(0 * u.deg, -45 * u.deg))
    assert not g.in_footprint(SkyCoord(0 * u.deg, 10 * u.deg))


class _FakeVizier(rx.get_catalog("nvss")):
    """NVSS with the network call replaced."""
    name = None  # don't register

    def load(self):
        pass

    def query(self, coord):
        from astropy.table import Table
        if coord.dec.deg > 50:
            return Table({"S1.4": [12.5]})[0]
        return None


def test_vizier_match_logic():
    m = _FakeVizier()
    assert m.match(src("det", 0, 60)) == {"NVSS_InFootprint": "Y", "NVSS_Flux_mJy": 12.5}
    assert m.match(src("nondet", 0, 0)) == {"NVSS_InFootprint": "Y", "NVSS_Flux_mJy": -999.0}
    r = m.match(src("south", 0, -60))
    assert r["NVSS_InFootprint"] == "N" and np.isnan(r["NVSS_Flux_mJy"])


def test_unknown_option_rejected(tmp_path):
    with pytest.raises(TypeError):
        rx.get_catalog("nvss")(bogus=1)


def test_config_end_to_end(tmp_path):
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "bobcat.csv").write_text("name,ra,dec\nA,10,10\n")
    (tmp_path / "data" / "list.txt").write_text("B\n")
    (tmp_path / "sources.csv").write_text("name,ra,dec,z\nA,10,10,0.1\nB,20,20,0.2\n")
    (tmp_path / "cfg.yaml").write_text(
        "match_radius_arcsec: 5\n"
        "input: {file: sources.csv}\n"
        "output: {file: out.csv}\n"
        "catalogs:\n"
        "  - {type: bobcat, file: data/bobcat.csv}\n"
        "  - {type: namelist, label: InList, file: data/list.txt}\n"
        "  - {type: nvss, enabled: false}\n"
    )
    df = rx.run_config(rx.load_config(str(tmp_path / "cfg.yaml")))
    assert list(df.columns) == ["Name", "z", "RA_deg", "Dec_deg", "BOBCat_Match", "InList"]
    back = pd.read_csv(tmp_path / "out.csv", comment="#")
    assert list(back["BOBCat_Match"]) == ["Y", "N"]
    assert list(back["InList"]) == ["N", "Y"]
    header = (tmp_path / "out.csv").read_text()
    assert "BOBCat_Match: Y/N: source is in BOBCat" in header


def test_column_clash_detected(tmp_path):
    (tmp_path / "l.txt").write_text("A\n")
    ms = [rx.get_catalog("namelist")(file=str(tmp_path / "l.txt"), label="X")] * 2
    with pytest.raises(ValueError):
        rx.crossmatch([src("A", 0, 0)], ms)
