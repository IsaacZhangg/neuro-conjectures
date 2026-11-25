# Catalogue

Every result with its exact scope. The README table is the short form.

## Model claims

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Hopfield network storage capacity](results/hopfield-capacity/) | Hebbian Hopfield network, ±1 neurons, unbiased random patterns, zero-temperature asynchronous dynamics started at a stored pattern; N = 500 to 4000; capacity defined as the load where the probability of ending with overlap ≥ 0.9 falls below one half; extrapolated to large N with an assumed N^-1/2 correction. | Simulation, fixed seed, recorded output. Replication of Amit, Gutfreund and Sompolinsky (1985). Not yet independently reviewed. |

## Theorems

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Rigorous Hopfield capacity bound](results/hopfield-rigorous-bound/) | Hebbian couplings. ±1 patterns: an energy barrier at Hamming radius ⌈0.0339N⌉ around every pattern, with a stable fixed point of overlap ≥ 0.932, probability 1 − e^{−cN}, for α < 0.12948 (Theorem 2, proof of the universality step written, review pending). Gaussian interference: the same barrier (Theorem 1), and a fixed point for α < 0.134716 (Theorem 3). Novelty audit: Zhou (2026, unrefereed) claims fixed points up to 0.1382, which supersedes Theorem 3; the barrier bound is new against the 0.08 on record. | Written proofs (Sudakov–Fernique comparison, capped Lindeberg, concentration); constants certified with exact rational arithmetic, cross-checked with Arb. Agent-drafted; ±1 proof under independent agent review. |

## Empirical claims

None yet. Screened targets are in [docs/candidates.md](docs/candidates.md).
