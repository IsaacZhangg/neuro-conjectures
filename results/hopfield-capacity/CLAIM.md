# Claim

## Statement tested

For the Hopfield model with Hebbian couplings J_ij = (1/N) Σ_μ ξ_i^μ ξ_j^μ
(i ≠ j), storing P = αN independent, unbiased random patterns ξ^μ ∈ {−1, +1}^N,
zero-temperature dynamics has stable states highly correlated with the stored
patterns only for α < α_c, with α_c ≈ 0.138 as N → ∞.

## Source of the value

Amit, Gutfreund and Sompolinsky (1985) derive α_c = 0.138 in the
replica-symmetric theory, with a jump of the error fraction from about 1.6% to
about 50% at the transition, and note that replica-symmetry breaking may raise
α_c slightly.

## Operational definitions

- Retrieval: starting exactly at a stored pattern, asynchronous dynamics ends
  with overlap m = (1/N) Σ_i s_i ξ_i^μ ≥ 0.9.
- Finite-size capacity α_c(N): load where the retrieval probability crosses
  1/2, by linear interpolation on the grid α = 0.10, 0.11, …, 0.24.
- Large-N capacity: intercept of a least-squares fit α_c(N) = α_c + c·N^-1/2.

## Decision rule (fixed before the full run)

Supported if the extrapolated α_c lies in [0.118, 0.158] (0.138 ± 0.02).
Refuted if it lies outside that interval and the result is stable to doubling
the number of samples.
