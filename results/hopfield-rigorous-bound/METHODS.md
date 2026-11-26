# Methods

## Model or data

Hebbian Hopfield couplings J_ij = (1/N) Σ_μ ξ^μ_i ξ^μ_j, i ≠ j, J_ii = 0.
Model G: one ±1 target pattern and M Gaussian interference patterns.
Model R (standard): all P patterns independent uniform ±1. Full definitions
in PROOF.md §1.

## Analysis plan

1. Reduce the energy barrier at Hamming radius k to an inequality between
   max_{σ∈Σ_k} ‖Gσ‖ and ‖G1‖ (PROOF.md §3).
2. Bound the expected maximum with the Sudakov–Fernique inequality
   (PROOF.md §4, Lemma 4).
3. Add Gaussian Lipschitz concentration and a χ² lower tail (PROOF.md §5).
4. Certify the constant: check.py computes rational intervals for π
   (Machin series), √(2π), exp(−τ²/2), Q(τ) (Taylor series of the error
   function) and A(τ) at τ₀ = 18259/10000, each series with an explicit
   geometric tail bound. The decision rule passes if the lower end of the
   A(τ₀) interval exceeds 0.113.

## Uncertainty

The certified interval has width below 10⁻¹⁴, from the √ precision
(40 digits) and the series tolerance (10⁻⁴⁵). No sampling is involved in the
decision.

## Robustness checks

- Independent recomputation of A(τ₀) with Arb ball arithmetic
  (python-flint, 200 bits), when installed. It agrees.
- Floating-point scan of sup_τ A(τ) (0.12948998 at τ ≈ 1.82587), and
  comparison with Stojnic's reported expression at x = τ*/√2 (agreement to
  about 3·10⁻⁹).
- Exploratory zero-temperature dynamics from the target, models G and R,
  N = 1000, α ∈ {0.10, 0.12, 0.129}. Final overlaps stay near 0.99. This
  is consistent with the theorem but proves nothing, and dynamics from the
  pattern is a weaker notion than the energy barrier.

## Deviations from plan

- Theorem 3 (reference-state refinement, 0.134716) was added after the
  decision rule. It was found while pushing the result further and is
  reported separately in the CLAIM amendment. The original rule and the
  Partial status are unchanged.

- The bound was derived before CLAIM.md was committed, so the decision rule
  was written after the constant was known, though before PROOF.md and
  check.py existed. The rule is the obvious one (beat 0.113) and was not
  tuned.
- The planned target was the ±1 model. Only the Gaussian-interference model
  was proved; the ±1 case is reduced to hypothesis (U). Hence Partial.
