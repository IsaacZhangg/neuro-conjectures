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

None yet. Screened targets are in [docs/candidates.md](docs/candidates.md).
