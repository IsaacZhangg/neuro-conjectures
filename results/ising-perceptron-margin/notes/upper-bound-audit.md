# Audit: Huang's upper-bound conditions at κ > 0 (agent audit, 2026-11-25)

Floating point unless marked [flint]. Written by a Claude agent; independent review pending.

- No hidden κ = 0 dependence outside Appendix B, where it is explicit (eqs (100), (109),
  Claims B.7–B.12). Typos that matter at κ ≠ 0: eq. (14) must read +α E log Ψ; the §7.5 F′
  identity omits κ (general form F′ = −F(F − (κ−x)/(1−q0))); Lemma 7.17's criterion should read
  α⋆κ E[N]. General-κ d0 = −(1−q0)ψ0 + α⋆κ E[N] replaces eq. (109).
- κ = 0 reproduces the paper: Condition 3.3 = 0.54436, λ0 = −0.19100, Fig. 1.

| κ | α⋆ | Cond 3.3 | λ0 | ∇²S̄⋆(1,0) eigenvalues |
|---|---|---|---|---|
| 0.25 | 0.6092595657 | 0.5175 | −0.2290 | −0.04559, −0.004458 |
| 0.5 | 0.4556645645 | 0.4895 | −0.2709 | −0.03077, −0.004345 |
| 0.75 | 0.3478935313 | 0.4611 | −0.3148 | −0.01857, −0.003962 |
| 1.0 | 0.2707000923 | 0.4333 | −0.3590 | −0.00972, −0.002692 |

- Condition 1.3 (global, two variables): grid scans of S̄⋆ (201² in tanh coordinates plus
  λ ∈ [−12, 12]²) show no positive value at κ = 0.25, 0.5, 0.75, 1.0; the unique maximum is
  (1, 0); best value at distance ≥ 0.1 is about −2.5e-5. At κ = 1 a second local maximum of S̄⋆
  appears at (1.141, 0.423) with S̄⋆ = −1.64e-6 (S⋆ = −0.0094 there).
- Failure points found above the tested range: S̄⋆ ≤ 0 fails at κ ≈ 1.0012 (S⋆ itself stays
  negative, so Condition 1.3 is not shown to fail); Claim 2.6 fails at κ ≈ 1.18.
- [flint] Claim 2.6 (∇²S̄⋆(1,0) negative definite) certified on boxes around the parameters at
  κ = 0.25, 0.5, 0.75, 1.0 (κ = 0.5: H11 = −0.0252 ± 4e-5, det = 1.34e-4 ± 7.5e-7), relying on
  the agent's transcription of the Lemma B.9 general-κ formula (finite differences agree to
  5e-7). The boxes must be checked to contain the certified parameter enclosures from
  code/certify_1d.py; the tail constant used there is the agent's own and needs review.
Code: upper-bound/huang_float/.
