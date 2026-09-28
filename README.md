# cosmic_infrastructure

## radio_crossmatch

Give it a list of extragalactic sources (names and RA/Dec) and it cross-matches every source against a set of radio and binary-supermassive-black-hole catalogs. You get back one CSV table with a header that describes every column.

Each catalog is its own small module, so adding a new survey means adding one file, and sometimes only a few lines of config.

### Built-in catalogs

| `type` | Catalog | Source | Output columns |
|---|---|---|---|
| `bigmac` | BIG MAC binary/dual AGN candidates | local CSV | `BIGMAC_detected`, `BIGMAC_{Primary,Secondary,Tertiary}Type` |
| `bobcat` | BOBCat binary SMBH candidates | local CSV (`name, ra, dec`) | `BOBCat_Match` |
| `rfc` | Radio Fundamental Catalog (VLBI) | local fixed-width file | `RFC_Match`, `RFC_F{s,l}{S,C,X,U,K}` (Jy) |
| `vlass` | VLASS epoch 1 (Gordon+2021) | VizieR `J/ApJS/255/30/comp` | `VLASS_Flux_mJy` |
| `first` | FIRST 14Dec17 | VizieR `VIII/92/first14` | `FIRST_InFootprint`, `FIRST_Flux_mJy` |
| `nvss` | NVSS | VizieR `VIII/65/nvss` | `NVSS_InFootprint`, `NVSS_Flux_mJy` |
| `vizier` | **any** VizieR table, set up entirely in the config | VizieR | `<prefix>_InFootprint`, `<prefix>_Flux_mJy`, extras |
| `namelist` | Y/N: is the name in a plain list? (e.g. the Nyland samples) | local text file | `<label>` |

To see what's installed, run `radio-crossmatch --list`.

### Install

```bash
git clone https://github.com/xiggystardust/cosmic_infrastructure.git
cd cosmic_infrastructure
pip install -e ".[test]"
```

### Quick start

```bash
radio-crossmatch examples/basic/config.yaml
```

A config lists the input, the output, and the catalogs to run, in the order you want their columns:

```yaml
match_radius_arcsec: 10          # default radius for every catalog
on_error: warn                   # "raise" stops on the first failed query
input:
  file: sources.csv              # columns: name, ra, dec (deg); extra columns are kept
output:
  file: crossmatched.csv
catalogs:
  - type: nvss
  - type: rfc
    file: data/rfc_2026b_cat.txt
  - type: vizier                 # a new survey, no code needed
    catalog: VIII/81B/sumss212
    prefix: SUMSS
    flux_column: St
    dec_max: -30
    radius_arcsec: 20            # per-catalog override
```

Relative paths resolve against `data_dir` if you set one, and otherwise against the folder the config is in. To switch a catalog off without deleting it, add `enabled: false` to its entry.

To read the output back in: `pd.read_csv("crossmatched.csv", comment="#")`.

### Python API

```python
import radio_crossmatch as rx

sources = rx.read_sources("sources.csv")
matchers = [rx.get_catalog("nvss")(), rx.get_catalog("rfc")(file="rfc_2026b_cat.txt")]
df = rx.crossmatch(sources, matchers)
rx.write_output(df, "out.csv", rx.column_descriptions(matchers))
```

### Adding a catalog

For a VizieR table, try the `vizier` type in your config first. For anything else, copy `src/radio_crossmatch/catalogs/_template.py` to a new file in the same folder and fill it in. If the class has a `name`, it registers itself and you can use it as `type: <name>`. [CONTRIBUTING.md](CONTRIBUTING.md) walks through the details.

### Matching conventions

- **Nearest match within the radius.** Local catalogs use the closest entry. VizieR queries ask for rows sorted by distance (`+_r`) and keep the first one.
- **Name matching** (BIG MAC, BOBCat) runs before position matching. Names are compared after upper-casing and removing whitespace. The default `name_match: contains` reproduces the original script, but it is loose: `NGC410` also matches `NGC4101`. Use `name_match: exact` to avoid that, or `none` to match on position only.
- **Footprints are approximate** declination cuts (FIRST: dec > −10°, NVSS/VLASS: dec > −40°). A non-detection inside the footprint is reported as `-999`. FIRST in particular does not cover the whole sky north of −10°, so treat its `-999` values with caution.
- **VizieR table IDs are checked** when the catalog loads. A wrong ID otherwise looks exactly like "no detection" for every source.

### Examples

- `examples/basic/`: a minimal config and a three-source input list.
- `examples/horlaville/`: the Horlaville et al. PTA-target analysis this package grew out of. It shows how to add a custom input format, exclusion lists and output splits on top of the generic pipeline.

### Tests

```bash
pytest -m "not network"   # offline
pytest                    # also runs live VizieR queries
```

### License

MIT; see [LICENSE](LICENSE).
