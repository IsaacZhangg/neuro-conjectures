# Catalogue

Every result with its exact scope. The README table is the short form.

## Model claims

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Hopfield network storage capacity](results/hopfield-capacity/) | Hebbian Hopfield network, ±1 neurons, unbiased random patterns, zero-temperature asynchronous dynamics started at a stored pattern; N = 500 to 4000; capacity defined as the load where the probability of ending with overlap ≥ 0.9 falls below one half; extrapolated to large N with an assumed N^-1/2 correction. | Simulation, fixed seed, recorded output. Replication of Amit, Gutfreund and Sompolinsky (1985). Not yet independently reviewed. |

## Theorems

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Rigorous Hopfield capacity bound](results/hopfield-rigorous-bound/) | Hebbian couplings, one ±1 target pattern, M Gaussian interference patterns: energy barrier at Hamming radius ⌈0.0339N⌉ and a stable fixed point with overlap ≥ 0.932, probability 1 − e^{−cN}, for α < 0.12948; a stable fixed point within distance ⌈0.0261N⌉ for α < 0.134716 (Theorem 3, reference-state refinement). The standard ±1 model is covered only conditionally on a universality hypothesis (U); the rigorous ±1 bound stays at 0.113. | Written proof (Sudakov–Fernique comparison and concentration); constant certified with exact rational arithmetic, cross-checked with Arb. Agent-drafted, not yet independently reviewed. |

## Empirical claims

None yet. Screened targets are in [docs/candidates.md](docs/candidates.md).
