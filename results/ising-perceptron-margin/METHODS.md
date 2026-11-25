# Methods (frozen plan, written before any certified numbers)

Claim under test: for margin κ > 0, the capacity of the ±1-weight perceptron
equals the Krauth–Mézard prediction α*(κ).

## Openness check (done 2026-11-25)

Searched the literature for a proof at κ ≠ 0. Found: Ding–Sun (arXiv 1809.07742,
lower bound, κ = 0), Huang (arXiv 2404.18902, upper bound; Theorem 3.6 holds for
every κ but only under numerical Conditions 1.3, 3.1, 3.3, 3.4 and Claim 2.6,
verified only at κ = 0), Shmalo (arXiv 2607.11010, certifies the κ = 0 conditions),
Nakajima (arXiv 2512.23195, uniqueness of the RS saddle for κ ≥ 0). No proof for
κ > 0 found, so the claim is treated as open.

## Plan

1. Upper bound at fixed κ > 0: Huang's Theorem 3.6 reduces it to numerical
   conditions. Certify each one-dimensional condition (3.1, 3.3, 3.4) and the
   local Hessian condition (Claim 2.6) with Arb ball arithmetic (python-flint).
   Condition 1.3 is a global two-dimensional inequality; certify it if feasible,
   otherwise give floating-point evidence and name it as a gap.
2. Lower bound at fixed κ > 0: audit every step of Ding–Sun that uses κ = 0,
   then certify (or give evidence for) Ding–Sun Condition 1.2 at κ.
3. κ values: 0 (sanity check against the literature), then 0.1, 0.25, 0.5, 0.75, 1.0.

## Decision rule (fixed in advance)

- Status Proved only if every condition in both directions is certified and every
  κ = 0-specific step in the cited proofs has a written replacement.
- Status Partial if some conditions are certified and the rest are named gaps.
- Never call a floating-point check a certification.

## Rigour of the numerics

Gaussian expectations are computed as Arb integrals over [−14, 14] (or [0, 14] for
even integrands) plus an explicit tail bound from a polynomial envelope of the
integrand. Envelopes used: ℰ(x) ≤ |x| + 1, 0 < ℰ′ < 1, |log Ψ(x)| ≤ x²/2 + |x| + 2,
log 2ch(y) ≤ |y| + log 2.
