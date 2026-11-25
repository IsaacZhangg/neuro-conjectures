"""Download the public data this result uses and verify checksums.

The data are not committed. They are taken from two published Python
packages on PyPI, which bundle the original files:

- wormneuroatlas 0.0.7.3 (Randi et al. 2023 signal-propagation atlas and the
  anatomical connectome it uses). Only data files are extracted; no package
  code is imported.
- cect 0.3.5 (OpenWorm C. elegans Connectome Toolbox), which bundles the
  Ripoll-Sanchez et al. 2023 neuropeptide connectome CSVs.

Usage: python fetch_data.py [--dest data]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
import zipfile
from pathlib import Path

WHEELS = {
    "wormneuroatlas": {
        "version": "0.0.7.3",
        "file": "wormneuroatlas-0.0.7.3-py3-none-any.whl",
        "sha256": "166289031784c9565f6b9205e30dc718bafbc092db3b46f1aef73d93dc4d1325",
        "extract": [
            "wormneuroatlas/data/funatlas.h5",
            "wormneuroatlas/data/aconnectome_default.h5",
            "wormneuroatlas/data/neuron_ids.txt",
            "wormneuroatlas/data/aconnectome_ids_ganglia.json",
            "wormneuroatlas/data/anatlas_neuron_positions.txt",
        ],
    },
    "cect": {
        "version": "0.3.5",
        "file": "cect-0.3.5-py3-none-any.whl",
        "sha256": "ea47f7a6b0d69a8d0bd3b0ade4037b7cc32d31ad0449f823d7ec153d6fd25057",
        "extract": [
            "cect/data/01022024_neuropeptide_connectome_short_range_model.csv",
            "cect/data/01022024_neuropeptide_connectome_mid_range_model.csv",
            "cect/data/01022024_neuropeptide_connectome_long_range_model.csv",
        ],
    },
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def wheel_url(name: str, version: str, filename: str) -> str:
    with urllib.request.urlopen(f"https://pypi.org/pypi/{name}/{version}/json", timeout=60) as r:
        meta = json.load(r)
    for u in meta["urls"]:
        if u["filename"] == filename:
            return u["url"]
    raise RuntimeError(f"{filename} not listed on PyPI")


def fetch(dest: Path) -> dict:
    dest.mkdir(parents=True, exist_ok=True)
    out = {}
    for name, w in WHEELS.items():
        whl = dest / w["file"]
        if not whl.exists() or sha256(whl) != w["sha256"]:
            url = wheel_url(name, w["version"], w["file"])
            print(f"downloading {w['file']}")
            urllib.request.urlretrieve(url, whl)
        got = sha256(whl)
        if got != w["sha256"]:
            raise RuntimeError(f"checksum mismatch for {w['file']}: {got}")
        with zipfile.ZipFile(whl) as z:
            for member in w["extract"]:
                target = dest / Path(member).name
                if not target.exists():
                    target.write_bytes(z.read(member))
                out[Path(member).name] = sha256(target)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default=str(Path(__file__).parent / "data"))
    args = ap.parse_args()
    for k, v in fetch(Path(args.dest)).items():
        print(k, v)
    return 0


if __name__ == "__main__":
    sys.exit(main())
