# Claim

## Statement tested

Cuntz, Mathy and Häusser (2012) derive that a tree which connects n target
points spread through a volume V with near-minimal cable has total length

    L ≈ (3/4π)^(1/3) · V^(1/3) · n^(2/3)        (their Eq. 1)

and, because branch points bp are proportional to n, that **L ∝ V^(1/3) · bp^(2/3)**.
They tested the branch-point form on light-microscopy reconstructions and
report a fitted exponent of 0.72 ± 0.10 across datasets (0.70 pooled).

This result tests the branch-point form, with exponent 2/3, on whole-cell
electron-microscopy reconstructions from two connectomes:

- **FlyWire** adult *Drosophila* whole brain (materialization 783), whole-neuron
  skeletons, one female fly.
- **MICrONS** mouse visual cortex (minnie65, skeleton release v661), dendritic
  skeletons of the cells typed in the Allen V1 column, one mouse.

The exact quantity, our reading of the paper's regression: the slope **β** of

    log( L / V^(1/3) )  =  c + β · log(bp)

fitted over the included neurons of one dataset. The claim predicts β = 2/3.

## Where it comes from

Cuntz H, Mathy A, Häusser M. A scaling law derived from optimal dendritic
wiring. PNAS 109(27):11014–11018, 2012. doi:10.1073/pnas.1200430109.
Eq. 1 (derivation in Results); branch-point form and the 0.72 ± 0.10 fit in the
NeuroMorpho analysis (Fig. 3 and text). Read on PMC (PMC3390826). The equation
itself is an image there that we could not load; the constant (3/4π)^(1/3) is
taken from the derivation text. See SOURCES.md.

## Operational definitions

- **Tree.** FlyWire: the public skeleton of the whole neuron (soma, neurite,
  axon and dendrite together, since FlyWire skeletons carry no compartment
  labels). MICrONS: SWC nodes labelled soma (1), basal dendrite (3) or apical
  dendrite (4); axon nodes (2) are dropped.
- **Pruning.** One pass removes every terminal segment (tip to the first branch
  point) shorter than τ = 1 µm, to suppress skeletonization twigs.
- **L.** Sum of edge lengths of the pruned tree, µm.
- **bp.** Number of vertices of degree ≥ 3 in the pruned tree.
- **V.** Convex hull volume of the pruned tree's vertices, µm³ (the paper's
  choice).

## Decision rule (frozen before the analysis reads the data)

Window W = [0.567, 0.767], that is 2/3 ± 0.1. The width covers the paper's own
spread (0.72 ± 0.10) and the finite-size value our estimator gives on
wiring-optimal benchmark trees (minimum spanning trees of random points,
checked in `check.py`, about 0.60).

For each dataset separately, with β and its 95% bootstrap interval from the
primary analysis in METHODS.md:

- **Supported in that dataset** if the whole interval lies inside W.
- **Refuted in that dataset** if the interval lies entirely outside W, and the
  interval also lies entirely outside W in each of the three robustness
  variants (τ = 0, τ = 5 µm, tips instead of branch points).
- **Contested in that dataset** if the primary analysis is decisive (Supported
  or Refuted) but a robustness variant reaches the opposite decisive verdict.
- **Inconclusive in that dataset** otherwise.

Status of the result (one word, set from the two dataset verdicts):

- Both Supported: **Partial** (two single-animal datasets; the paper's claim is
  general, so this is a narrower scope than stated).
- Both Refuted: **Refuted** for these two datasets and this form of the law.
- Any dataset Contested: **Contested**.
- Any other combination: **Partial**, with each dataset's verdict stated in
  the README.

The check fails (exit 1) without setting a status if the estimator misses W
on the minimum-spanning-tree benchmark, because the test would then not be
able to recognise a wiring-optimal tree.
