# Follow-up: does physical proximity explain the pair-specific peptide gain?

Added 2026-11-25, after the main run, and committed before any of the
analysis below was run. It does not change METHODS.md, CLAIM.md or the
status those set. It asks one question about the main result: does the small
pair-specific part of the no-wired-path gain (G_none = +0.111 against a
degree-preserving null mean of +0.081, so about +0.03) survive controls for
how close the two neurons are?

Why: the only peptide model that beat its null, the Ripoll-Sánchez
short-range model, links neurons whose processes share a nerve bundle. Nearby
neurons may respond to each other for reasons that have nothing to do with
peptides (shared bundles, cross-talk, unmapped contacts). If so, the +0.03
would vanish once proximity is held fixed.

What was looked at before this file was frozen: the format of the two
proximity files (a header line of neuron names and one x, y, z line per
neuron; a JSON map from ganglion name to neurons). No proximity value was
compared with any response or peptide link.

## Proximity data

- Soma positions: `anatlas_neuron_positions.txt` from wormneuroatlas 0.0.7.3
  (first line: neuron names; line k + 1: x, y, z of the k-th name). Distance
  d(i, j) is Euclidean. Units do not matter, because only quantile bins are
  used. Names are standardised as for the peptide CSV (AWCL → AWCOFF,
  AWCR → AWCON). Neurons without a position get no distance.
- Ganglion: `aconnectome_ids_ganglia.json` from the same package. Every key
  other than "head" and "pharynx" names a ganglion and lists its neurons.
  same_ganglion(i, j) is true if both neurons are listed under one ganglion.
- Proximity class of a pair: distance quintile (quintiles over all ordered
  neuron pairs i ≠ j that both have positions) × same_ganglion, so 10 classes.
  Pairs where either neuron lacks a position form class −1.

## Fixed settings

Primary network (short-range), standard labels (q < 0.05), two-step score
S(A) = A + A·A, no-wired-path split as in the main analysis, seed 20261125,
2000 cluster-bootstrap resamples over stimulated neurons, 1000 nulls per test.

## F1: proximity as a covariate

AUROC is computed within proximity classes only (the same stratified AUROC
as R4: sum of Mann–Whitney U over classes divided by the sum of
n_responding × n_not_responding; class −1 kept as its own class). This removes
any effect of proximity on responding, as binned. Compute G_none and G_all
this way for the real peptide layer and for 1000 ordinary degree-preserving
rewirings (the main analysis's null), also scored within classes.

## F2: shuffle only among comparable distances

A degree-preserving rewiring in which each double-edge swap
(r1, c1), (r2, c2) → (r1, c2), (r2, c1) is accepted only if
class(r1, c2) = class(r1, c1) and class(r2, c1) = class(r2, c2). Each neuron
keeps its in- and out-degree and its number of outgoing-to-incoming links in
each proximity class per row, and the network keeps how many links fall in
each class. 10 × (number of links) accepted swaps, at most
1000 × (number of links) attempts; the fraction of links moved is recorded.
Score with the ordinary (unstratified) AUROC, as in the main analysis.

## Decision for this follow-up

For each of F1 and F2, the pair-specific gain survives if the real G_none
exceeds the null with one-sided p = (1 + number of null values ≥ real) /
1001 below 0.05. Report the excess (real G_none minus null mean) with a 95%
cluster-bootstrap interval of the real value minus that fixed null mean.

- Survives proximity: survives both F1 and F2.
- Mixed: survives exactly one.
- Explained by proximity: survives neither.

G_all is reported the same way but is not part of this decision.

## Exploratory, not part of the decision

- AUROC of −d (nearer = higher score) for responding pairs, all pairs and
  no-wired-path pairs.
- Share of short-range peptide links in each distance quintile.

## Deviations from this follow-up plan

None yet.
