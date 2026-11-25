# Claim

## Target, frozen word for word

From the candidate list (entry 30, "All-pattern capacity of dense associative
memory with Gaussian patterns"), the line that decides whether this counts as
solved:

> **Counts as solved if:** There is a proof that α_c equals the union bound (the large-λ case α_c = (log 2)/4 alone is a clean first result), or a proof that α_c is strictly larger than the union bound.

## Statement tested

Let ξ^1, …, ξ^P be independent standard Gaussian vectors in R^N (entries iid
N(0,1), so E‖ξ‖² = N), with P = ⌊e^{αN}⌋ and α > 0 fixed as N → ∞.

For λ ∈ (0, ∞) and each pattern μ define the noise term

    F_λ^μ = (1/λ) log Σ_{ν ≠ μ} exp(λ ξ^μ·ξ^ν),

and for λ = ∞ let F_∞^μ = max_{ν ≠ μ} ξ^μ·ξ^ν. Pattern μ is **retrieved** when
the signal beats the noise at x = ξ^μ:

    ‖ξ^μ‖² > F_λ^μ.

Let R_N(α, λ) be the event that every pattern μ = 1, …, P is retrieved, and

    α_c(λ) = sup{ α : P(R_N(α, λ)) → 1 as N → ∞ }.

Let α_c^lb(λ) be the union-bound threshold of Lucibello and Mézard, eq. (10):
α_c^lb(λ) = sup{α : A(α, λ) − α > 0}, with A from their eq. (9) and the
Gaussian rate functions of their Appendix D (eqs. 63 to 67).

**Claim (target):** α_c(λ) = α_c^lb(λ). In the large-λ case this reads
α_c = (log 2)/4.

We aim to prove the sharp form: for α < α_c^lb(λ), P(R_N) → 1, and for
α > α_c^lb(λ), P(R_N) → 0.

## Where it comes from

Lucibello, C. and Mézard, M. "Exponential Capacity of Dense Associative
Memories." Phys. Rev. Lett. 132, 077301 (2024); arXiv:2304.14964.

- Energy and dynamics: eqs. (1), (2). Signal versus noise split: eqs. (3), (4), (6).
- Union bound and the definition of α_c^lb: Sec. III, eqs. (8) to (10).
- Gaussian rate functions: Appendix D, eqs. (63) to (67); "for Gaussian patterns
  there is a gap at any λ".
- Large-λ criterion: Appendix F.2, eq. (86), ‖ξ^μ‖² > ξ^μ·ξ^ν for all ν ≠ μ.
  The paper: "It remains an open question whether we have α_c^lb = α_c or not."

## Operational definitions (our reading)

The paper defines retrieval dynamically (gradient iterations from a nearby start
converge near the pattern), but its union bound, and so α_c^lb, is computed for
the static signal-dominance event ‖ξ^μ‖² > F_λ^μ of eq. (6), and its large-λ
numerics use eq. (86), which is the λ = ∞ case of the same event. We take that
static event as the definition of "retrieved". This is our reading, stated here
before any proof work. The dynamical version is not covered by this claim and
any link to it will be stated separately as a remark.

## Decision rule (committed before the proof work)

- **Proved** if PROOF.md gives a complete argument that for every λ ∈ (0, ∞],
  P(R_N) → 1 for α < α_c^lb(λ) and P(R_N) → 0 for α > α_c^lb(λ), every numerical
  constant used in it is certified (exact arithmetic or interval arithmetic),
  and an independent reviewer (labelled as an agent if it is one) finds no gap.
- **Partial** if the argument covers only λ = ∞, or only a stated range of λ,
  or only one side.
- **Disproved** (of the equality) if the argument shows α_c(λ) > α_c^lb(λ)
  for some λ.
- Numerical brute-force runs at small N are supporting evidence only and cannot
  set the status.
