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

None yet. Screened targets are in [docs/candidates.md](docs/candidates.md).
