# Methods

Pre-registered plan. This file and CLAIM.md were committed before any code read
NeuroMorpho.org metadata or reconstructions. Field names below come from the
NeuroMorpho.org API reference page (read 2026-11-25), not from the data.

## Data

- Source: NeuroMorpho.org REST API, `https://neuromorpho.org/api/neuron`,
  all pages at `size=500`. The download date and the API's `totalElements`
  are recorded in `check-output.json`; NeuroMorpho has no versioned snapshot,
  so the neuron id list and per-file SHA-256 checksums of what was analysed
  are recorded instead (`manifest.json`, committed after the run; no raw files
  are committed).
- Files: the standardised reconstructions,
  `https://neuromorpho.org/dableFiles/<archive, lower case>/CNG%20version/<neuron_name>.CNG.swc`.
  (This path follows common usage of the site; it is not in the API reference.)
- License: NeuroMorpho.Org content is CC BY 4.0. Each analysed neuron's archive
  and `reference_pmid` are recorded so the original papers can be cited.

## Inclusion, by metadata (before download)

A neuron is eligible if all hold:

1. `attributes` starts with `Diameter` (excludes `No Diameter`).
2. `domain` contains `Dendrites` and does not contain `No Dendrites`.
3. `experiment_condition` is `Control`.
4. `cell_type` contains no entry with `glia` (case-insensitive).

Groups are (archive, species, cell_type list as given). A group is eligible if
it has at least 10 eligible neurons. From each eligible group, up to 20 neurons
are drawn uniformly without replacement with NumPy seed 20261125 (groups in
sorted order, neurons sorted by `neuron_id` before drawing).

## Inclusion, by file (after download, before any exponent is computed)

A downloaded neuron is analysed only if its dendritic samples (types 3 and 4):

1. number at least 50;
2. take at least 10 distinct diameter values;
3. have a modal diameter covering at most 50% of dendritic samples (screens
   out default or constant diameters);
4. give at least 10 valid bifurcations (definition in CLAIM.md; all three
   diameters > 0).

A group qualifies for the decision rule if at least 10 of its drawn neurons
pass. Neurons that fail to download, fail to parse, or fail a filter are
counted by reason in the output.

## Analysis plan

1. Parse each CNG SWC file; build the tree from parent ids.
2. For each valid bifurcation, compute d₀, d₁, d₂ with the 5 µm windows in
   CLAIM.md.
3. Per cell, find η̂ by evaluating the loss on a grid η = 0.500, 0.505, …,
   6.000 and taking the minimum (no SciPy dependency). Record how many cells
   sit at a grid bound.
4. Per group, median η̂ and its 95% bootstrap interval over cells
   (2000 resamples).
5. Pooled estimate and interval as in CLAIM.md.
6. Apply the decision rule.

Bootstrap generator: NumPy `default_rng(20261125)`, used in a fixed order.

## Uncertainty

The group interval resamples cells. Cells from the same animal are not
independent, and NeuroMorpho metadata has no animal identifier, so group
intervals are likely too narrow. The pooled interval resamples archives (labs),
which is the coarsest independent unit available.

## Robustness checks (each re-runs steps 2 to 6)

- R1: section-mean diameters (mean over every sample of the parent and daughter
  sections, as in Desai-Chowdhry et al. 2022) instead of 5 µm windows.
- R2: 10 µm windows.
- R3: drop bifurcations where d₀ = d₁ = d₂ exactly (likely quantisation or
  default values).

If any of R1 to R3 gives a different status from the primary analysis, the
status is Contested.

## Estimator validation (synthetic, run in quick mode)

Synthetic binary trees with d₀^1.5 = d₁^1.5 + d₂^1.5 exactly, random daughter
splits, 10% lognormal multiplicative noise on every sample diameter and 0.1 µm
quantisation. The estimator must return a median η̂ within ±0.10 of 1.5 over
20 synthetic cells, and within ±0.10 of 2.0 for trees built with η = 2. This
checks the code, not the claim; the observed bias is reported.

## Secondary analyses (reported, do not set the status)

- Proximal against distal: per group, η̂ fitted separately on bifurcations of
  centrifugal order 1 and 2 against order 3 and higher (Desmond and Levy 1984
  report a proximal/distal difference in dentate granule cells; we could only
  read part of their abstract).
- Fraction of bifurcations per cell where d₀^(3/2) is within ±10% of
  d₁^(3/2) + d₂^(3/2).

## Amendment A1 (2026-11-25, committed before any Blue Brain file is read)

This adds a secondary test. It does not change the primary analysis, its
data, or how the status is set. It exists because NeuroMorpho.org is not
reachable from the cloud machine this repo's agent runs on, while the Blue
Brain Open Data bucket is. Only the bucket's file listing (names and sizes)
had been seen when this section was written; no morphology file had been
opened.

Data:

- Source: Blue Brain Open Data, public AWS S3 bucket `openbluebrain`
  (us-west-2), documented at https://github.com/BlueBrain/OpenData and in the
  AWS Open Data Registry entry `bluebrain_opendata`. License CC BY 4.0
  ("Copyright (c) 2024 Blue Brain Project/EPFL licensed under CC-BY 4.0").
- Files: every `.swc` file under
  `Experimental_Data/Reconstructed_morphologies/Categorized/Neurons/Rat/Cortex/SSCX/<m-type>/downloaded_from_nexus/`.
  The listing on 2026-11-25 showed 1074 such files in 60 m-type folders, 40
  of them with at least 10 files. These are rat somatosensory cortex
  reconstructions; we take the folder name as the cell-type label.
- Other Blue Brain folders (mouse, human, `Categorized_clean`, `.asc` and
  `.h5` files) are not used. They mix sources and formats, and SWC is the only
  format `check.py` parses.
- A file name that appears in more than one m-type folder is dropped from all
  of them, so that no cell is counted twice.
- Files are cached under `data/bluebrain/` (not committed). The S3 key, size
  and SHA-256 of each analysed file go in `manifest-bluebrain.json`.

Analysis: as the primary plan, steps 1 to 4, with the same file filters
(Inclusion, by file), windows, grid, robustness checks R1 to R3 and secondary
analyses. Up to 20 files per m-type are drawn with the same seed, from keys
sorted by name. A group is one m-type and qualifies with at least 10 passing
cells.

Pooled estimate: median of group medians, with a 95% percentile bootstrap
over groups (2000 resamples). All files come from one institute, so there is
no between-lab resampling; this interval is likely too narrow.

Reporting: the decision rule's labels are applied and the resulting word is
reported as the "Blue Brain secondary label". It never sets or changes the
result's status. If it disagrees with the primary status, the README says so
next to the primary result. Output goes to `check-output-bluebrain.json`.

Known limits of A1: one institute, one species, one brain area; the files are
as Blue Brain published them through its Nexus store, and we have not checked
whether their diameters were edited by a repair or curation step.

## Known limits before the run

- Diameters near the light-microscope resolution limit (≈0.3 to 0.5 µm) are
  quantised and noisy, which can bias η̂ upwards (daughters look as thick as
  parents). R3 and the secondary distal split partly address this.
- Per-bifurcation exact solutions for η are ill-conditioned (Liao, Liang and
  Howard 2021 report a mean of 7.2 and SD of 8.8 that way), which is why the
  per-cell least-squares fit on log diameters is used instead.
- NeuroMorpho's cell-type labels and the `attributes` flag are as submitted;
  we do not re-check them.

## Deviations from plan

- 2026-11-26, code fix, no change to the rule: the first full run (on Isaac's
  Mac, 2026-11-25) matched no neuron, because the NeuroMorpho API returns
  `experiment_condition` as a list (`['Control']` for 370,683 of 495,348
  neurons) and `check.py` compared it with the string `Control`. The check now
  keeps a neuron only when every listed condition is `Control`, which is the
  inclusion rule above. Found from a count of field values only; no
  reconstruction file had been downloaded.

Otherwise none to the primary plan. Amendment A1 (above) adds a secondary test on a
second dataset, committed before that dataset was read.
