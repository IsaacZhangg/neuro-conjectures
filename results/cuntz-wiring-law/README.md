# Cuntz wiring law on EM connectomes

Does total dendrite length scale as V^(1/3) · bp^(2/3), the wiring-optimal
scaling law of Cuntz, Mathy and Häusser (2012), in whole-cell EM
reconstructions from the fly and mouse connectomes?

Status: Contested
Kind: Empirical claim
Reproducibility: 🟡 runnable check with recorded output, independent AI agent review only

## Result

The fitted exponent is about 0.53 in both connectomes, below the predicted
2/3 and below the window 2/3 ± 0.1 that the frozen decision rule allows.

| Dataset | Neurons | β, primary (τ = 1 µm) | No pruning | 5 µm pruning | Tips as n | Verdict |
| :-- | --: | :-- | :-- | :-- | :-- | :-- |
| FlyWire 783, whole neurons | 6,025 (402 types) | 0.537 [0.532, 0.541] | 0.556 [0.551, 0.561] | 0.459 [0.454, 0.464] | 0.561 [0.556, 0.565] | Refuted |
| MICrONS v661, dendrites | 1,309 (E and I) | 0.528 [0.514, 0.541] | 0.523 [0.510, 0.537] | 0.533 [0.521, 0.545] | 0.613 [0.598, 0.628] | Contested |

β is the slope of log(L / V^(1/3)) on log(bp); brackets are 95% bootstrap
intervals over neurons within cell types. The window is [0.567, 0.767].

- **FlyWire:** every interval, in the primary analysis and all three
  robustness variants, lies below the window, so the rule gives Refuted.
  The margin is thin: the tips variant ends 0.001 below the window edge, and
  an equally valid draw of the 15 neurons per type gives Inconclusive (see
  METHODS.md, deviation 3). Read this as Refuted, marginally.
- **MICrONS:** the primary interval lies below the window, but counting tips
  instead of branch points puts it inside, so the rule gives Contested.
- **Status:** the frozen rule gives Contested when any dataset is Contested.

The estimator gives 0.597 on wiring-optimal benchmark trees (minimum spanning
trees of random points), so part of the gap is the finite-size bias of the
method. The data sit below even that benchmark.

Secondary analyses (they do not set the status):

- Free fit, log L = c + a·log V + b·log bp: FlyWire a = 0.201, b = 0.604;
  MICrONS a = 0.279, b = 0.539. The claim predicts a = 1/3, b = 2/3. Both
  datasets fall short on both exponents.
- Per cell type, FlyWire (15 neurons each): β = 0.50 ± 0.17 (mean ± SD over
  402 types, median 0.49); 28% of types fall inside the window. Cuntz and
  colleagues report 0.72 ± 0.10 over NeuroMorpho cell types.
- MICrONS by class: excitatory 0.530, inhibitory 0.515.

Full numbers are in [check-output.json](check-output.json); per-neuron
measurements in [neuron-stats.csv](neuron-stats.csv).

## Files

- [CLAIM.md](CLAIM.md): exact statement and decision rule.
- [SOURCES.md](SOURCES.md): primary source and later work.
- [METHODS.md](METHODS.md): data, analysis plan, robustness checks, deviations.
- [check.py](check.py) and [wiring.py](wiring.py): entry point and tree
  measurements. `python check.py --quick` (benchmark, no network) /
  `python check.py` (full run, about 8 minutes; downloads skeletons).
- [check-output.json](check-output.json): recorded full run.
- [manifest.json](manifest.json): sources and SHA-256 of every skeleton used.
- [neuron-stats.csv](neuron-stats.csv): derived per-neuron L, bp, tips and V.
  FlyWire rows are CC BY-NC 4.0, MICrONS rows CC BY 4.0; see
  [DATA-LICENSE.md](DATA-LICENSE.md).
- [REVIEW.md](REVIEW.md): review record.

## Scope and limits

- One fly (FlyWire 783) and one mouse (MICrONS minnie65 v661). The bootstrap
  is over neurons within each animal, not over animals.
- This tests the branch-point form with a fitted exponent. The synapse form
  with fixed exponents was tested on EM data by Ferreira Castro and Cardona
  (2024), and is not retested here.
- FlyWire skeletons mix axon and dendrite, and many fly neurons spread over
  several neuropils, so the convex hull overstates the volume the arbor fills.
  FlyWire reconstructions also miss some fine twigs.
- The result depends on how branch points are counted: pruning length moves
  the FlyWire β from 0.46 to 0.56, and using tips moves the MICrONS β to 0.61.
  This sensitivity is why MICrONS is Contested.
- MICrONS skeletons are dendrites only (axons removed in the release) and are
  fairly coarse (median 21 branch points), which may undercount branching.
- A low exponent here does not show that neurons waste wire. It shows that,
  measured this way on EM skeletons, length grows more slowly with branch
  points than the 2/3 law predicts.
