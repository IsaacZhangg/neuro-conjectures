# Audit: does the Ding–Sun lower bound carry over to κ > 0? (agent audit, 2026-11-25)

Floating point only. Not a proof. Written by a Claude agent; independent review pending.

Source: Ding and Sun, arXiv 1809.07742 v1 (Sections 1, 2, 7, 8 read in full; 3–6 and 9–10
read for statements and every use of κ), plus targeted parts of the journal version
(PTRF 193 (2025) 627–715), which drops Sections 8–10.

## Verdict

No genuine obstacle found. What must change for κ > 0:

- Prop 1.3 / §7 computer-assisted constants (eq. (153), Lemmas 7.4–7.5, Cors 7.6–7.7):
  recompute the parameter boxes at each κ.
- Lemma 7.3: its proof is wrong as written for κ ≠ 0. The reflection z̄ = 2κ√q − z does not
  preserve φ, and there is a sign slip (ξ(z̄) − ξ(z) = −2qζ). It can be repaired for κ ≥ 0
  because φ(z̄) ≤ φ(z) on {z ≤ κ√q}. Not repaired this way for κ < 0.
- §§8–10 (the analytic control near λ = 0 and λ = 1 and the bulk sweep) are written with
  "we fix κ ≡ 0": ξ = γz, odd symmetry in z, splits at z = 0. They need rewriting with
  ξ = γz + κ/√(1−q), splitting on the sign of ξ; monotonicity in q flips at z = κ√q.
- The CLT step of Theorem 1.4 uses (H+P)″(0) < 0, which is stronger than S″(0) < 0; state it as
  an explicit assumption and check it.
- The κ_all → κ rescaling needs joint continuity of q⋆ in (κ, α) (implicit function theorem).

Everything else (Lemmas 2.4, 2.9, 7.1, Cor 7.2, Sections 3–6, the identities S(0) = S′(0) = 0,
S(1) = −G⋆, S′ → ∞ as λ ↑ 1) is written for general κ.

## Condition 1.2 (journal form)

S⋆(λ) < 0 for all λ ∉ {0, 1} in [λ_min, 1], with S⋆″(0) < 0 and lim_{λ↑1} S⋆′(λ) > 0.

## Floating-point evidence at α = α⋆(κ)

| κ | α⋆ | λ_min | max S on \|λ\| ≥ 0.1 | (H+P)″(0) | S″(0) | S(λ_min) |
|---|---|---|---|---|---|---|
| 0 | 0.8330786 | −0.4243 | −1.30e-3 | −0.0596 | −0.2771 | −0.0523 |
| 0.1 | 0.7330167 | −0.4214 | −1.28e-3 | −0.0625 | −0.2744 | −0.0504 |
| 0.25 | 0.6092596 | −0.4174 | −1.28e-3 | −0.0676 | −0.2704 | −0.0478 |
| 0.5 | 0.4556646 | −0.4118 | −1.25e-3 | −0.0767 | −0.2633 | −0.0440 |
| 1.0 | 0.2707001 | −0.4036 | −1.17e-3 | −0.0938 | −0.2464 | −0.0376 |

S ≤ −0.066 λ² on |λ| ≤ 0.9 at every κ tested; S(1−ι) ≈ −c√ι with c ≈ 0.42–0.52, so S′ → +∞.
For κ ≤ 0.25 the s-tilt term A is needed on a negative-λ window; for κ ≥ 0.5 H + P < 0 alone.
Code: code/dingsun_float/ (ds_S.py, analyze.py; run `python analyze.py <κ> 2001`).
