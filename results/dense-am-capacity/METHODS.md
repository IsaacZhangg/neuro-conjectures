# Methods

## Model

As in CLAIM.md: P = ⌊e^{αN}⌋ iid N(0, I_N) patterns, interaction strength
λ ∈ (0, ∞], static retrieval event ‖ξ^μ‖² > F_λ^μ for every μ.

## Proof plan (frozen before proof work)

Notation: s_μ = ‖ξ^μ‖²/N; I(s) = (s − 1 − log s)/2 is the large-deviation rate
of s_μ (Lucibello and Mézard eq. 63 with s = r²).

Two facts drive everything:

1. A Gaussian vector is its norm times an independent uniform direction, so
   for any fixed unit vector u, ξ^ν·u ~ N(0,1) exactly, whatever ‖ξ^ν‖ is.
2. Failure of "all retrieved" is monotone: adding patterns can only add
   competitors, so P(R_N) is non-increasing in P.

**Lower bound (α < α_c^lb ⇒ P(R_N) → 1).** Union bound over μ, then condition on
s_μ. Given ξ^μ, the overlaps ξ^μ·ξ^ν/N are iid N(0, s_μ/N). Bound the upper tail
of F_λ^μ by splitting the competitors into overlap levels on a fine grid: levels
whose expected count is exponentially large have counts concentrated
(Chernoff); levels with vanishing expected count contribute at most a few
competitors except with super-exponentially small probability. This gives
F_λ^μ/N ≤ max(φ_typ(s_μ), largest overlap/N) + o(1) outside a negligible event,
where φ_typ is the typical REM free energy (Lucibello and Mézard eq. 5). The
union bound then reproduces A(α, λ) − α > 0.

**Upper bound (α > α_c^lb ⇒ P(R_N) → 0).** Split the patterns into halves A
(candidates) and B (competitors). Two failure mechanisms:

- (a) "short pattern swamped by typical noise": some candidate has s_μ below the
  root s_0(α, λ) of φ_typ(s) = s. Choose such a candidate by its norm only; its
  noise from B is bounded below by one overlap level whose count concentrates.
- (b) "one competitor beats a short pattern": condition on A. The events
  "competitor ν ∈ B beats some short candidate" are independent across ν given
  A, so P(no failure | A) ≤ exp(−|B| w(A)), where w(A) is the probability that
  one fresh Gaussian vector beats some short candidate. Lower-bound w(A) by the
  first two Bonferroni terms: the first is a sum over candidates of a function of
  their norms only (concentrates); the second is controlled in expectation
  through the bivariate Gaussian tail and the density of the angle between two
  independent uniform directions. At s = 1/2 this gives failure for
  2α > (log 2)/2, i.e. α > (log 2)/4, using monotonicity to stay in α < 1/4
  where the Bonferroni term is small.

Then check that the two mechanisms together are exactly the complement of
A(α, λ) − α > 0.

## Constants

- The large-λ threshold (log 2)/4 is closed form (minimum of I(s) + s/2 at
  s = 1/2).
- The λ-crossover below which mechanism (a) binds, and any other number used in
  an inequality, are certified with ball arithmetic (python-flint / Arb) or
  exact rationals.
- Lean 4: formalize small calculus lemmas only if it is cheap in this
  environment (Lean and Mathlib are not preinstalled here).

## Numerical support (does not set the status)

Brute force for λ = ∞ at small N: sample P = ⌊e^{αN}⌋ Gaussian patterns, check
eq. (86) for all pairs, estimate P(R_N) over a grid of α near (log 2)/4 and
several N with a fixed seed. Also compare the λ-dependence of α_c^lb computed
from the paper's formula with the threshold predicted by mechanisms (a) and (b).

## Uncertainty

Not applicable to the proof. Brute-force probabilities get Wilson intervals over
independent seeds.

## Robustness checks

Recompute every constant independently (a separate agent review recomputes the
exponents and the crossover). Check the λ = ∞ statement against the paper's own
Fig. 8 numerics.

## Deviations from plan

- Lean: not done. No Lean toolchain here and the needed probability (Cramér
  bounds, Gaussian tails, angle densities) is thin in Mathlib, so it was not
  cheap. Constants are certified by Arb interval arithmetic instead.
- Fig. 8 comparison: not done. The brute-force check at λ = ∞ (N ≤ 40) is
  reported in check-output.json as support only.
- The proof turned out not to need a separate shell-and-δ argument for the
  lower bound beyond Lemma 1 and a cell covering; the plan's steps are
  otherwise as frozen.
