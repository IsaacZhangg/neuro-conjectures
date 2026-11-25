# Catalogue

Every result with its exact scope. The README table is the short form.

## Model claims

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Hopfield network storage capacity](results/hopfield-capacity/) | Hebbian Hopfield network, ±1 neurons, unbiased random patterns, zero-temperature asynchronous dynamics started at a stored pattern; N = 500 to 4000; capacity defined as the load where the probability of ending with overlap ≥ 0.9 falls below one half; extrapolated to large N with an assumed N^-1/2 correction. | Simulation, fixed seed, recorded output. Replication of Amit, Gutfreund and Sompolinsky (1985). Not yet independently reviewed. |

## Theorems

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [All-pattern capacity of dense associative memory, Gaussian patterns](results/dense-am-capacity/) | Dense associative memory with log-sum-exp energy at interaction strength λ ∈ (0, ∞], P = e^{αN} iid standard Gaussian patterns, static signal-dominance retrieval (Lucibello and Mézard eqs. 6 and 86), N → ∞. Target: all-pattern threshold equals their union bound, (log 2)/4 at large λ. | Proved: α_c(λ) = min(α_a(λ), (log 2)/4) for every λ ∈ (0, ∞], sharp transition; (log 2)/4 for λ ≥ 0.70092. Agent-reviewed. |

## Empirical claims

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Neuropeptide links and C. elegans signal propagation](results/celegans-neuropeptides/) | Wild-type Randi et al. 2023 signal-propagation atlas (wormneuroatlas 0.0.7.3; 23,902 measured pairs, q < 0.05); wired connectome bundled with it; Ripoll-Sánchez et al. 2023 peptide connectome (short-range primary, mid and long as checks); binary path score A + A·A scored by AUROC; cluster bootstrap over stimulated neurons; 1000 degree-preserving peptide rewirings. | Pre-registered reanalysis, recorded output. Partial: worse over all pairs, better (AUROC 0.61) for pairs with no wired path within two steps, mostly explained by peptide degree; mid- and long-range models add nothing beyond degree. A frozen follow-up finds the small pair-specific part is mostly explained by soma proximity (Mixed). Agent review only. |

Screened targets are in [docs/candidates.md](docs/candidates.md).
