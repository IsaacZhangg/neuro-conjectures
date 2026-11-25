# Catalogue

Every result with its exact scope. The README table is the short form.

## Model claims

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Hopfield network storage capacity](results/hopfield-capacity/) | Hebbian Hopfield network, ±1 neurons, unbiased random patterns, zero-temperature asynchronous dynamics started at a stored pattern; N = 500 to 4000; capacity defined as the load where the probability of ending with overlap ≥ 0.9 falls below one half; extrapolated to large N with an assumed N^-1/2 correction. | Simulation, fixed seed, recorded output. Replication of Amit, Gutfreund and Sompolinsky (1985). Not yet independently reviewed. |

## Theorems

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Tree committee machine capacity, K = 3](results/tree-committee-k3/) | Three sign hidden units on disjoint thirds of the input, majority output; i.i.d. Gaussian patterns, labels fixed or independent; storage of all ⌊αn⌋ patterns with probability → 1. Proves α_c(3) ≥ 1.55 (open gap 1.55 to 3.43). Lower bound only; ±1 inputs not covered. | Written proof (three-stage construction, convex Gaussian min-max theorem, Cover's theorem) and Arb-certified constants with recorded output. Agent-drafted, not yet independently reviewed. |

## Empirical claims

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Neuropeptide links and C. elegans signal propagation](results/celegans-neuropeptides/) | Wild-type Randi et al. 2023 signal-propagation atlas (wormneuroatlas 0.0.7.3; 23,902 measured pairs, q < 0.05); wired connectome bundled with it; Ripoll-Sánchez et al. 2023 peptide connectome (short-range primary, mid and long as checks); binary path score A + A·A scored by AUROC; cluster bootstrap over stimulated neurons; 1000 degree-preserving peptide rewirings. | Pre-registered reanalysis, recorded output. Partial: worse over all pairs, better (AUROC 0.61) for pairs with no wired path within two steps, mostly explained by peptide degree; mid- and long-range models add nothing beyond degree. A frozen follow-up finds the small pair-specific part is mostly explained by soma proximity (Mixed). Agent review only. |

Screened targets are in [docs/candidates.md](docs/candidates.md).
