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

None yet.
