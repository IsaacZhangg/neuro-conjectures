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
| [SCS random network global stability](results/scs-global-stability/) | dx/dt = −x + J tanh(x), J_ij iid N(0, g²/N). Proved: origin globally asymptotically stable for all initial conditions with probability ≥ 1 − exp(−N(1 − √2 g)²/(2g²)) when g < 1/√2. Open: 1/√2 ≤ g < 1. | Written proof (Lyapunov function plus Gaussian random-matrix bounds), numerical spot checks with recorded output. Independent agent review. Novelty: partly known; same method and constant published for random Lotka–Volterra (see NOVELTY.md). |
| [Rigorous Hopfield capacity bound](results/hopfield-rigorous-bound/) | Hebbian couplings. ±1 patterns: an energy barrier at Hamming radius ⌈0.0339N⌉ around every pattern, with a stable fixed point of overlap ≥ 0.932, probability 1 − e^{−cN}, for α < 0.12948 (Theorem 2, proof of the universality step written, passed an independent agent review). Gaussian interference: the same barrier (Theorem 1), and a fixed point for α < 0.134716 (Theorem 3). Novelty audit: Zhou (2026, unrefereed) claims fixed points up to 0.1382, which supersedes Theorem 3; the barrier bound is new against the 0.08 on record. The barrier (energy above the pattern's own on the whole sphere) is a different, stronger notion than Zhou's AGS capacity (an energy well near the pattern), so the two numbers are not comparable bounds on one quantity. | Written proofs (Sudakov–Fernique comparison, capped Lindeberg, concentration); constants certified with exact rational arithmetic, cross-checked with Arb. Agent-drafted; ±1 proof passed an independent agent review (2026-11-26); not yet checked by a person. |
| [Tree committee machine capacity, K = 3](results/tree-committee-k3/) | Three sign hidden units on disjoint thirds of the input, majority output; i.i.d. Gaussian patterns, labels fixed or independent; storage of all ⌊αn⌋ patterns with probability → 1. Proves α_c(3) ≥ 1.55 (open gap 1.55 to 3.43). Lower bound only; ±1 inputs not covered. | Written proof (three-stage construction, convex Gaussian min-max theorem, Cover's theorem) and Arb-certified constants with recorded output. Agent-drafted, not yet independently reviewed. |

## Empirical claims

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Neuropeptide links and C. elegans signal propagation](results/celegans-neuropeptides/) | Wild-type Randi et al. 2023 signal-propagation atlas (wormneuroatlas 0.0.7.3; 23,902 measured pairs, q < 0.05); wired connectome bundled with it; Ripoll-Sánchez et al. 2023 peptide connectome (short-range primary, mid and long as checks); binary path score A + A·A scored by AUROC; cluster bootstrap over stimulated neurons; 1000 degree-preserving peptide rewirings. | Pre-registered reanalysis, recorded output. Partial: worse over all pairs, better (AUROC 0.61) for pairs with no wired path within two steps, mostly explained by peptide degree; mid- and long-range models add nothing beyond degree. A frozen follow-up finds the small pair-specific part is mostly explained by soma proximity (Mixed). Novelty: partly known (wired baseline and peptide negative in Sale & Richardson 2026). Agent review only. |

Screened targets are in [docs/candidates.md](docs/candidates.md).
