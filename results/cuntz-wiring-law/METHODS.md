# Methods

Plan committed before any code computed L, bp or V from the connectome data.
What was opened before this commit, to learn file formats: bucket listings,
the cell-type and annotation tables (to count candidate neurons per type), one
MICrONS SWC file and its metadata JSON, and the vertex array of one FlyWire
skeleton (shape and units only). No length, branch count or volume was
computed from any real neuron before the commit.

## Data

### FlyWire (fly)

- Skeletons: public precomputed skeletons `skeletons_mip_1` in the Google
  Cloud bucket `flywire_v141_m783` (FlyWire materialization 783), read with
  `cloud-volume` over HTTPS. Units nm, converted to µm.
- Annotations: `flyconnectome/flywire_annotations`,
  `supplemental_files/Supplemental_file1_neuron_annotations.tsv`, commit
  a83b2776d60d5764cef36b927f5f9679c16c47a2.
- Inclusion: `super_class` in {central, optic, visual_projection,
  visual_centrifugal, endocrine}; `status` empty (drops `outlier_seg` and
  `outlier_bio`); `cell_type` given. Sensory, motor, ascending and descending
  neurons are left out because they leave the brain and are cut at its edge.
- Sampling: cell types with at least 15 eligible neurons (407 types). From
  each, 15 neurons drawn without replacement with
  `numpy.random.default_rng(20261125)`, types in sorted order (6,105 neurons).

### MICrONS (mouse)

- Skeletons: SWC files in the BossDB open-data bucket,
  `s3://bossdb-open-data/iarpa_microns/minnie/minnie65/skeletons/v661/skeletons/`,
  named `<root_id>_<nucleus_id>.swc`, units µm. Read over HTTPS.
- Cell types: `iarpa_microns/minnie/minnie65/cell_types/allen_visp_column_soma_coarse_types_v1.csv`
  (nucleus id, valid flag, classification system, class). Included: rows with
  class `excitatory` or `inhibitory`, matched to a skeleton by nucleus id. If a
  nucleus id matches more than one skeleton file, the neuron is dropped.
- Strata for resampling: excitatory, inhibitory.

### Both

- After pruning, a neuron is excluded if L < 100 µm, bp < 10 or V ≤ 0, or if
  its skeleton cannot be read. Counts of every exclusion are recorded.
- No raw data is committed. `check.py` downloads what it needs. It writes
  `neuron-stats.csv` (one row per included neuron: dataset, id, stratum, L,
  bp, tips, V for each pruning setting) and `manifest.json` (source paths,
  versions, SHA-256 of each downloaded file or skeleton array).

## Analysis plan

1. For every included neuron and τ in {0, 1, 5} µm, prune (one pass), then
   measure L, bp, tips and V as defined in CLAIM.md.
2. Primary analysis per dataset, τ = 1 µm: ordinary least squares of
   y = log(L / V^(1/3)) on x = log(bp) over all included neurons. β is the
   slope.
3. 95% interval: percentile bootstrap, 2,000 resamples, resampling neurons
   with replacement within strata (FlyWire: cell type; MICrONS: excitatory or
   inhibitory), seed 20261125.
4. Apply the decision rule in CLAIM.md.

## Benchmark (estimator calibration)

Minimum spanning trees of uniform random points in boxes of random aspect
ratio, with 50 to 3,000 points and sides 50 to 500 µm (150 trees, seed
20261125). The same pipeline (τ = 0) gives β on these trees. MSTs are
wiring-optimal trees without a root constraint, so their β is what the
estimator reports for trees that follow the law at these finite sizes.

## Uncertainty

Each dataset is one animal, so the bootstrap is over neurons within one
animal. The interval describes variation across neurons of that brain, not
across animals.

## Robustness checks

Three variants of the primary analysis, each with its own interval:

- R1: τ = 0 (no pruning).
- R2: τ = 5 µm.
- R3: tips instead of branch points as n (τ = 1 µm).

How they enter the verdict is set in CLAIM.md.

## Secondary analyses (reported, do not set the status)

- Free-exponent fit per dataset: log L = c + a·log V + b·log bp, with
  intervals (same bootstrap). The claim predicts a = 1/3, b = 2/3.
- Per cell type (FlyWire, 15 neurons each) and per class (MICrONS): β by type,
  summarised as mean ± SD across types, to compare with the paper's
  0.72 ± 0.10 across NeuroMorpho datasets.
- The constant: median of L / (V^(1/3) · bp^(2/3)) per dataset.

## Known limits planned for

- FlyWire skeletons mix axon and dendrite in one tree, and fly neurons often
  spread across several neuropils, so the convex hull includes empty space.
  The claim's V is the volume the arbor spans; the hull overstates it for such
  cells.
- FlyWire reconstructions miss some fine twigs (Ferreira Castro and Cardona
  estimate about 47% dendritic completeness by synapse attachment), which lowers
  bp and L.
- MICrONS cells near the volume edge are cut; the V1 column sits in the
  middle of the volume, which limits this.
- n is branch points, not synapses. The synapse form is the one Ferreira
  Castro and Cardona tested; this result does not retest it.

## Deviations from plan

1. **MICrONS matching (before any MICrONS neuron was measured).** The first
   field of `allen_visp_column_soma_coarse_types_v1.csv` turned out to be an
   annotation id, not a nucleus id, so no row matched a skeleton by id. Rows
   are now matched to nucleus ids through
   `iarpa_microns/minnie/minnie65/nucleus_detection/nucleus_detection_v0.csv`
   by exact equality of the soma point (x, y, z in voxels); a row with no
   unique exact match is dropped and counted. All 1,347 E/I rows matched
   exactly, and 1,342 of those nuclei have a v661 skeleton. Found during a
   smoke run on five neurons per dataset; the five FlyWire neurons measured in
   that run were not looked at beyond checking that the code ran.
2. **FlyWire count (correction of a number, not of the rule).** The plan
   said 407 cell types and 6,105 neurons. That count was made before the
   `status` filter. With the filter, as the plan specifies, 402 types have at
   least 15 eligible neurons, so the run used 402 × 15 = 6,030 neurons. The
   sampling rule itself is unchanged.

The following points were found by the independent review after the full run.
None changes a verdict; they are recorded so the plan and the code can be
compared exactly.

3. **Order before sampling.** The plan did not say how a type's root ids are
   ordered before `rng.choice`. `check.py` sorts them ascending. With file
   order instead, the FlyWire sample shares only 2,207 of 6,030 neurons; that
   draw gives primary β = 0.538 [0.534, 0.543] and the tips variant
   0.564 [0.559, 0.568], whose upper end crosses the window edge, so FlyWire
   would be Inconclusive rather than Refuted. The overall status is Contested
   either way.
4. **Exclusions are applied once.** The size rule (L ≥ 100 µm, bp ≥ 10,
   V > 0) is applied at τ = 1 µm, and the same neurons enter every variant.
   The plan did not say this. Applying the rule separately at each τ moves
   FlyWire R2 from 0.459 to 0.470; no verdict changes.
5. **Benchmark recipe.** The code draws the number of points and the box side
   log-uniformly, and the two shorter sides are 0.3 to 1 times the long side,
   so they can be as short as 15 µm, below the stated 50 to 500 µm. The
   written recipe, implemented by the reviewer, gives β = 0.609 instead of
   0.597; both are inside the window.
