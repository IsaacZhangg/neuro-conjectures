# Licenses of the derived data in this folder

The code here is MIT, like the rest of the repository. The derived data files
carry the terms of the datasets they come from.

| File | What it holds | Terms |
| :-- | :-- | :-- |
| neuron-stats.csv, rows with `dataset = MICrONS` | Per-neuron length, branch points, tips and hull volume computed from MICrONS minnie65 v661 skeletons | CC BY 4.0. Cite the MICrONS Consortium (Nature 2025). |
| neuron-stats.csv, rows with `dataset = FlyWire` | The same measurements and the cell-type label, computed from FlyWire 783 skeletons and annotations | **CC BY-NC 4.0** (non-commercial). Cite Dorkenwald et al. 2024 and Schlegel et al. 2024, plus the annotation papers listed in SOURCES.md. |
| manifest.json | Source paths, neuron ids and SHA-256 checksums | Ids and checksums only; same attribution as above. |
| check-output.json | Fitted slopes, intervals and summary statistics | Results of this check (repository text license), with attribution to both datasets. |

No raw skeletons or source tables are committed; `check.py` downloads them.
