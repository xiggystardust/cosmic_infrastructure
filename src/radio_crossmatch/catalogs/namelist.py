"""Flag sources whose name appears in a plain list of names.

Useful for "has this galaxy already been observed by programme X?" checks
(e.g. the Nyland et al. ATLAS3D and MASSIVE radio samples).  Config::

    - type: namelist
      label: NylandATLAS              # output column name
      file: data/nyland+2016-atlas3D.txt
      format: whitespace              # name = first whitespace token
    - type: namelist
      label: NylandMASSIVE
      file: data/eden_H26_crossmatch_sorted.txt
      format: fixed                   # name = characters [start:end]
      start: 0
      end: 13

Names are compared after upper-casing and removing whitespace.
"""

from ..catalog import CatalogMatcher
from ..utils import normalize_name, yn


class NameList(CatalogMatcher):
    name = "namelist"
    description = "Y/N flag for presence in a plain list of names (no positions)"

    def __init__(
        self,
        file: str,
        label: str,
        format: str = "whitespace",
        column: int = 0,
        start: int = 0,
        end: int = 13,
        description: str = "",
        **kwargs,
    ):
        super().__init__(**kwargs)
        if format not in ("whitespace", "fixed", "csv"):
            raise ValueError("format must be 'whitespace', 'fixed' or 'csv'")
        self.file, self.label, self.format = file, label, format
        self.column, self.start, self.end = column, start, end
        self.columns = {label: description or f"Y/N: name is listed in {file}"}
        self.names = set()

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, value):
        self._label = value

    def _extract(self, line):
        if self.format == "fixed":
            return line[self.start:self.end]
        parts = line.split(",") if self.format == "csv" else line.split()
        return parts[self.column] if len(parts) > self.column else ""

    def load(self):
        with open(self.file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                n = normalize_name(self._extract(line))
                if n:
                    self.names.add(n)

    def empty_result(self):
        return {self.label: "N"}

    def match(self, source):
        self.ensure_loaded()
        return {self.label: yn(normalize_name(source.name) in self.names)}
