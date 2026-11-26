# Ding–Sun second-moment exponent: computer-assisted bound on the bulk

**Agent-written (Claude), not independently reviewed.** Treat every claim here as
unverified until a person has checked the argument and the code.

## What is certified

For the Ising perceptron with margin κ, let S(λ) = H(λ) + P(λ) + A(λ) be the
Ding–Sun second-moment exponent (arXiv 1809.07742 v1, eqs (34)–(48)), with
κ_all = κ and (α, q, ψ) the replica-symmetric point at α = α⋆(κ).
The files `certificate_kappa<κ>*.json` give a finite list of boxes, with a
rigorous upper bound on S on each box, showing

    S(λ) < 0   for all λ in [λ_min, −δ] ∪ [δ, 1 − ε]

for the (δ, ε) recorded in each certificate, simultaneously for every
(α, q, ψ) in that certificate's parameter box (`params` field; `params_<κ>*_tight.json`).
Each parameter box contains (α⋆, q⋆(α⋆), ψ⋆(α⋆)) and (α, q⋆(α), ψ⋆(α)) for every α
in an interval around α⋆ (width 2·10⁻⁹ for the `_tight` boxes, 2·10⁻¹² for the
`_narrow_tight` boxes; see "Parameters").

## Results

| κ | δ | ε | parameter box | boxes | worst certified bound on S | where (λ) | wall time |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.05 | 0.01 | `params_0.05_tight.json` | 24 | −2.89·10⁻⁴ | [0.05, 0.068] | 5.4 min |
| 0.05 | 0.005 | 0.001 | `params_0.05_tight.json` | 41 (1 bisected) | −2.72·10⁻⁶ | [0.005, 0.0068] | 9.2 min |
| 0.13 | 0.005 | 0.001 | `params_0.13_tight.json` | 42 (2 bisected) | −2.78·10⁻⁶ | [0.005, 0.0068] | 12.4 min |
| 0.05 | 0.001 | 0.0001 | `params_0.05_narrow_tight.json` | 55 | −1.22·10⁻⁷ | [0.001, 0.00136] | 9.3 min |
| 0.13 | 0.001 | 0.0001 | `params_0.13_narrow_tight.json` | 55 | −1.21·10⁻⁷ | [0.001, 0.00136] | 9.5 min |

Every box of every certificate passed (no unresolved boxes). "Bisected" boxes
first failed (design too wide near λ ≈ 0.91–0.99) and both halves then passed;
the failed parents are listed under `failed_then_bisected` and are not part of
the cover. On the negative side the worst bound is at the box ending at −δ
(e.g. −1.30·10⁻⁷ at κ = 0.05, δ = 0.001). Near λ = 1 the last box (ending at 1 − ε)
has bound ≈ −4.7·10⁻³ (ε = 10⁻⁴), so ε is not the binding constraint. Wall times
are on 2 worker processes.

The certificates for δ ∈ {0.05, 0.005} were produced by an earlier revision of
`certify_bulk.py` that differs from the committed one only in two float design
choices (box half-width cap 0.6(1−λ) instead of 0.3(1−λ); λ-sub-ball width floor
2·10⁻⁴ instead of 10⁻⁷). The rigorous bound per box is the same code.

The neighbourhoods of λ = 0 and λ = 1 are not covered here; they are handled by
a separate analytic argument (another work item).

## Method

All arithmetic is Arb ball arithmetic (python-flint 0.9, 96-bit precision);
all integrals are `acb.integral` (rigorous Gauss–Legendre with analytic error
bounds) plus explicit tail bounds. Floating point (`ds_S.py`) is used only to
choose boxes and the parameter s; it plays no role in the bounds.

1. **Drop the infimum.** A(λ) = inf_s B(λ, s) ≤ B(λ, s) for any s. With s fixed
   per box, the I_0 terms of P and B cancel:

       S(λ) ≤ H(λ) + PB(λ),  PB(λ) = −P⋆ + ψ(1−q) f/2 + s²/2 − √ψ√(1−q) s √f + I_s(λ),
       f = (1−λ)/(1+λ).

2. **Parametrize by A.** λ = ℓ(A) = E D_{√ψZ}(A)/(1−q), with
   D_H(A) = (A²−1)c²/(Δ+1)² = c(Δ−1)/(Δ+1), c = 1−m², Δ = √(1+(A²−1)c), m = tanh H.
   Δ is increasing in A, so D and ℓ are increasing in A. Boxes are intervals
   [A0, A1] with shared endpoints, so their λ-images tile an interval; the first
   negative box starts at A = 0 (λ = λ_min = −E(1−|m|)²/(1−q)).

3. **H is concave in λ, with an explicit derivative.** For the pair law
   p = ¼((1+m)²+D, c−D, c−D, (1−m)²+D) one checks p11/p12 = (Δ+m)/(1−m) and
   p22/p12 = (Δ−m)/(1+m), so p11 p22 / p12² = (Δ²−m²)/c = A². Since
   ∂Γ/∂D = −¼ ln(p11 p22/p12²) = −½ ln A, we get dH/dA = −½ ln A · (1−q) ℓ′(A), so

       H′(λ) = −(1−q) ln A(λ) / 2,

   which is decreasing in λ. Hence H(λ) ≤ H(λc) + H′(λc)(λ − λc) on the whole box
   (λc = ℓ(Ac), Ac an interior point). Γ is evaluated from the weights
   (Δ+m)/(1−m), 1, 1, A²(1−m)/(Δ+m), which avoids cancellation at small and large A.

4. **Second-order bound for PB.** With an exact interior point λ⋆,
   PB(λ) ≤ PB(λ⋆) + PB′(λ⋆) d + ½ max(0, sup PB″) d², d = λ − λ⋆. PB(λ⋆) and
   PB′(λ⋆) are nested `acb.integral`s at a point (item 5). sup PB″ over the box
   only needs to be crude, so it is a **range (interval-Riemann) enclosure** in
   real ball arithmetic over λ-sub-balls of width ≤ min(0.02, 0.1(1 − λ)) (item 7).
   The final box bound is C + max(g d0 + M d0²/2, g d1 + M d1²/2) with
   C = H(Ac) + H′(λc)(λ⋆ − λc) + PB(λ⋆), g = H′(λc) + PB′(λ⋆) (convex in d, so the
   maximum is at an endpoint).

5. **The 2-D integral I_s.** Outer z ∈ [−10, 10] in unit pieces; inner
   ν ∈ [ξ_z, c] with c = ⌈max(sup ξ, 0)⌉ + 10 on each piece, written as
   ν = ξ + (c − ξ)t, t ∈ [0, 1], with weight φ(ν)/Ψ(ξ) = exp(−(ν²−ξ²)/2)/(√(2π) Ψ(ξ)e^{ξ²/2}).
   PB and PB′ are computed in one complex integral of g0 + i·g1 (both real on the
   real axis). Truncations are bounded by envelopes |g| ≤ e0 + e1|ν| + e2ν² using
   |ln Ψ(x)| ≤ x²/2 + |x| + 2, 0 < E(x) ≤ max(x,0) + 0.8, 0 < E′ < 1,
   E_{φ_ξ}|ν| ≤ |ξ| + 1.6, E_{φ_ξ}ν² = 1 + ξE(ξ) ≤ 1 + ξ² + 0.8|ξ|, and
   |ξ_z| ≤ κ/√(1−q) + |z|√q/√(1−q) (κ ≥ 0); see `dsarb.py::_envelope` and `Igen`.
   For the outer tail, the envelope is a polynomial in |z| of degree ≤ 2 with
   nonnegative coefficients, so it is ≤ P(L) z²/L² for |z| ≥ L.

6. **Numerics that matter for rigor.** Arb's erfc is very pessimistic on
   complex balls with large positive real part, so the scaled tail
   g(x) = Ψ(x)e^{x²/2} is enclosed by its value at the exact midpoint plus
   rad · (2π)^{−1/2}/c², using g(x) = (2π)^{−1/2} ∫₀^∞ e^{−xt−t²/2} dt and
   Re x ≥ c > ½ on the ball. Every function with a branch cut (log, sqrt) is
   called with the integrator's `analytic` flag, and inside the nested integral
   the inner integrand checks analyticity whenever either level requests it.

7. **Boundary-layer coordinates and the PB″ enclosure.** As λ → 1 the inner
   integrand in ν has a layer of width w = √(1−λ²) at ν = ξ and its λ-derivatives
   cancel catastrophically. With ν = ξ + w y (y ≥ 0),
   E_{ν∼φ_ξ} ln Ψ(arg) = E(ξ) ∫₀^∞ w e^Q h(x) dy, Q = −ξwy − w²y²/2,
   x = ξ r − λ y − c, r = √((1−λ)/(1+λ)), c = E(ξ)s/kk, h = ln Ψ. The y-limits do not
   depend on λ, so derivatives go under the integral with no cancellation
   (`Ip01_y` for PB, PB′ at λ⋆ ≥ 0.6; `I2_y` for PB″ on λ ≥ 0.6). For λ < 0.6 PB″
   uses ν = ξ + u (`I2_real`). The range enclosures use: on a z-cell the exact
   Gaussian mass Ψ(za) − Ψ(zb) times the integrand's range (mean value theorem,
   φ > 0); on each u- or y-cell the cell length times the range; E, E′ and ln Ψ on a
   real ball from monotonicity only (E increasing, E − x decreasing and positive,
   i.e. 0 < E′ < 1 (Sampford 1953); ln Ψ decreasing; Ψe^{x²/2} decreasing), with
   E′ = E·(E − x) enclosed as a product of the two monotone ranges. The y-tail
   y > Y uses |h|, |h′|, |h″| ≤ 2(t+1)φ(t) for x = −t ≤ −1, a polynomial majorant
   of the coefficients, and a geometric (ratio < 1, checked in Arb) bound.
   The z-tail |z| > 6 of PB″ uses the same ν-form envelope as item 5.

## Parameters

`params_<κ>.json` is the output of `certify_1d.certify` (Huang Conditions 3.1,
3.3, 3.4; the RS fixed point and α⋆ by sign changes), run with κ passed as the
decimal string, so `arb(κ)` encloses the exact decimal. `tighten_params.py`
replaces ψ0 = R(α-box, Q0) by the tighter box [R(α_lb, q⋆(α_lb)), R(α_ub, q⋆(α_ub))],
valid because q⋆(α) is increasing (T = P∘R increasing in α, sup ∂T/∂q < 1 on
the certified box) and ψ⋆ = P⁻¹(q⋆) with P increasing.

`params_<κ>_narrow.json` is the same certification run on an α-box of width
2·10⁻¹² centred at the root of the linear interpolation of G between the endpoints of
the wide box (and q-box ±10⁻¹⁰); it lies inside the wide α-box. The narrow box shrinks the
parameter-induced width of H(λ) from ~6·10⁻⁷ to ~5·10⁻¹⁰, which is what allows
δ = 0.001 (S(0.001) ≈ −1.4·10⁻⁷).

On whether α⋆ alone suffices: Ding–Sun's proof of Theorem 1.4 (Section 7 of v1,
the paragraph after the display that re-derives (47)) needs S_{0,α}(λ) uniquely maximized at 0 for α in a left
neighbourhood (α′, α⋆), and obtains it from Condition 1.2 at α⋆ plus
P″(0) + H″(0) < 0 and continuity. The certificate here holds on the whole
α-box, which includes α⋆ and a small two-sided neighbourhood, so it serves
either reading.

## Files

- `dsarb.py` — rigorous evaluation of ℓ, H, H′, P⋆, H⋆, λ_min; I_s and I_s′ at a point
  (`Igen('p01')` for λ⋆ < 0.6, `Ip01_y` for λ⋆ ≥ 0.6); PB″ on a λ-ball (`PB2_range`, via
  `I2_real` / `I2_y`). `Igen(2, ...)` and `I2_range` are earlier attempts at PB″ that
  the certificates do not use.
- `certify_bulk.py` — box design (float), per-box bound (Arb), bisection of failed boxes, coverage check.
- `tighten_params.py`, `run_cert.py` — parameter enclosures (`run_cert.py` imports
  `certify_1d.py` from `/home/claude/work`, which is not in this folder).
- `kmf.py` — float α⋆, q⋆, ψ⋆ (starting values only).
- `float_design/ds_S.py` — float reference, used only to place boxes and choose s.
- `params_*.json` — parameter certificates. `certificate_kappa<κ>_delta<δ>_eps<ε>.json` —
  every box with all enclosures and its bound, the coverage check, timings.
- `logs/` — run logs. `logs/old*_superseded.txt` are logs of aborted earlier runs with superseded PB″ code (one box there ended with bound +inf); they certify nothing.

## Reproduce

    export DS_FLOAT_DIR=$PWD/float_design
    /home/claude/venv-ising/bin/python run_cert.py 0.05 0.781074328272 0.57205355
    /home/claude/venv-ising/bin/python tighten_params.py params_0.05.json
    /home/claude/venv-ising/bin/python certify_bulk.py 0.05 params_0.05_tight.json 0.005 0.001 out d005e001
    # narrow box (centre from the wide box's G values) and the delta = 0.001 run
    /home/claude/venv-ising/bin/python run_cert.py 0.05 0.7810743282532177 0.572053536013 1e-12 1e-10 _narrow
    /home/claude/venv-ising/bin/python tighten_params.py params_0.05_narrow.json
    /home/claude/venv-ising/bin/python certify_bulk.py 0.05 params_0.05_narrow_tight.json 0.001 0.0001 out narrow
    # kappa = 0.13: run_cert.py 0.13 0.7059325315484992 0.5847018890457119 (wide),
    #               run_cert.py 0.13 0.7059325315101112 0.5847018482133 1e-12 1e-10 _narrow
