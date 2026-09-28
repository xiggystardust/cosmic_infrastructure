# Contributing a catalog module

This page covers how to add a new catalog to `radio_crossmatch`.

## 1. Pick the right starting point

| Your catalog is… | Do this |
|---|---|
| a VizieR table with one flux column | **No code.** Use `type: vizier` in your config (see `examples/basic/config.yaml`). |
| a VizieR table that needs custom logic or a default setup | Subclass `VizierCatalog` in `catalogs/vizier.py` (see `NVSS`, `FIRST`). |
| a file on disk with RA/Dec (and maybe names) | Subclass `LocalPositionalCatalog` (see `bobcat.py`, `bigmac.py`, `rfc.py`). |
| anything else (TAP, a web API, image cutouts, …) | Subclass `CatalogMatcher` and implement `match()`. |

## 2. Write the module

Copy `src/radio_crossmatch/catalogs/_template.py` to `catalogs/<yourcatalog>.py`. Every module in `catalogs/` whose name doesn't start with `_` is imported automatically.

A matcher needs:

```python
class MySurvey(CatalogMatcher):
    name = "mysurvey"            # registry key -> `type: mysurvey` in configs
    description = "One line, shown by `radio-crossmatch --list`"
    columns = {                  # ordered output columns + header descriptions
        "MYS_Match": "Y/N: detected in MySurvey",
        "MYS_Flux_mJy": "MySurvey 1.4 GHz flux density (mJy)",
    }

    def __init__(self, some_option=1, **kwargs):   # options come from the config entry
        super().__init__(**kwargs)                 # handles `radius`
        self.some_option = some_option

    def load(self):              # optional: one-off setup, called before matching
        ...

    def empty_result(self):      # optional: the row for "no match" (default: all NaN)
        return {"MYS_Match": "N", "MYS_Flux_mJy": float("nan")}

    def match(self, source):     # source.name, source.coord (SkyCoord or None), source.extra
        ...
        return {"MYS_Match": "Y", "MYS_Flux_mJy": 3.2}
```

Rules:

- `match()` must return exactly the keys in `columns`. The pipeline checks this.
- **Prefix your columns** with the catalog name. If two matchers produce the same column, the pipeline raises an error.
- Handle `source.coord is None` (no position available) by returning `empty_result()`.
- Use `self.radius` (an astropy Quantity) rather than a hard-coded radius.
- Raise exceptions on real failures instead of calling `exit()`. The user chooses `on_error: raise` or `warn`.
- Log with `logging.getLogger(__name__)` rather than `print`.
- If you need to override `match_all()` (for example to send one bulk query for the whole sample), keep its `on_error` behaviour.

Classes outside this package can register themselves too: subclass `CatalogMatcher` with a `name` and import your module before you build the config.

## 3. Test it

Add tests to `tests/`. Tests that don't touch the network should use small synthetic files written to `tmp_path` (see `tests/test_crossmatch.py`). Mark tests that hit a live service with `@pytest.mark.network`.

```bash
pytest -m "not network"
```

## 4. Document it

Add a row to the catalog table in `README.md`. Say where users can download the catalog file (if it's local), what its columns mean, and which units they're in.
