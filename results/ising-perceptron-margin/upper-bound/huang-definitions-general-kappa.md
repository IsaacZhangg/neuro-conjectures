# General-κ definitions in Huang, arXiv:2404.18902 (v2)

These were transcribed from the arXiv HTML (v2) and checked against my own derivations.
Equation numbers are the paper's. Notes marked **[note]** are my own; they are not proved.

Throughout, Z ~ N(0,1), th = tanh and ch = cosh.

## 1. Special functions (eqs (12), (13))

- **(12)** φ(x) = e^{-x²/2}/√(2π) and Ψ(x) = ∫_x^∞ φ.
- **(13)** E(x) = φ(x)/Ψ(x), and
  F_{1-q}(x) = (1-q)^{-1/2} E( (κ-x)/(1-q)^{1/2} ).
- Lemma 4.21 (valid for all real x):
  - 0 ≤ E(x) ≤ |x|+1
  - E' = E(E-x) ∈ (0,1)
  - E'' ∈ (0,1)
  - E''' ∈ (-1/2, 13)
- Consequence, for all x and κ:
  F'_{1-q}(x) = -(1-q)^{-1} E'(u), with u = (κ-x)/√(1-q).
  So -1/(1-q) < F' < 0.
- **[note]** In §7.5 the paper prints "F'_{1-q0}(x) = -F(F - x/(1-q0))". The general-κ form is
  F' = -F(F - (κ-x)/(1-q0)).
  The display just before it uses the κ form correctly.

## 2. Replica-symmetric threshold

- P(ψ) = E th(ψ^{1/2}Z)².
- R_α(q) = α E F_{1-q}(q^{1/2}Z)².
- **(14)** G(α,q,ψ) = -(1-q)ψ/2 + E log(2ch(ψ^{1/2}Z)) **+** α E log Ψ( (κ - q^{1/2}Z)/(1-q)^{1/2} ).
  - **[note]** (14) is printed with "−α E log Ψ". The proof of Lemma 7.15 identifies S⋆(Λ_{1,0}, √(1-q0)) with G(α⋆,q0,ψ0) using "+".
  - Condition 3.1 also needs G⋆ to be decreasing, and that holds only with "+", since log Ψ < 0.
  - I use "+". With it, α⋆(0) = 0.8330786 is reproduced.
- **Condition 3.1.** There exist 0<α_lb<α_ub and 0<q_lb<q_ub<1 such that, for every α in (α_lb,α_ub):
  - sup_{q∈(q_lb,q_ub)} (P∘R_α)'(q) < 1;
  - there is a unique q⋆(α,κ) in (q_lb,q_ub) with q⋆ = P(R_α(q⋆)); set ψ⋆ = R_α(q⋆);
  - G⋆(α) = G(α,q⋆,ψ⋆) is strictly decreasing, with a unique root α⋆(κ).
- Set q0 = q⋆(α⋆,κ), ψ0 = ψ⋆(α⋆,κ) and γ0 = q0/(1-q0).
- Prop 3.2 (from DS18, κ=0 only):
  - α ∈ (0.833078599, 0.833078600)
  - q ∈ (0.56394907949, 0.56394908030)
  - ψ ∈ [2.5763513100, 2.5763513224] (from App. B.1)

## 3. First-moment functional (§2.6, eqs (6)–(8))

- Random variables:
  - Ḣ ~ N(0,ψ0) and Ĥ ~ N(0,q0), on separate spaces; each expectation involves only one of them.
  - M = th(Ḣ).
  - N = F_{1-q0}(Ĥ) > 0.
- **(6)** ent(Λ) = E H((1+Λ(Ḣ))/2), where H(p) = -p log p - (1-p) log(1-p).
  - For Λ = th(h(Ḣ)): ent = E[log(2ch h) - h th h].
- **(7)** For s ≥ 0:
  S⋆(Λ,s) = ½ s² ψ0 + ent(Λ) + α⋆ E log Ψ{ U_Λ + sN },
  where
  U_Λ = (κ - (E[MΛ]/q0) Ĥ - (E[ḢΛ]/ψ0) N) / √(1 - E[MΛ]²/q0).
- **(8)** S⋆(Λ) = inf_{s≥0} S⋆(Λ,s).
  - Λ_{λ1,λ2}(x) = th(λ1 x + λ2 th x).
  - S⋆(λ1,λ2) = S⋆(Λ_{λ1,λ2}).
- **S̄⋆(λ1,λ2) = S⋆(Λ_{λ1,λ2}, √(1-q0))** (§2.6, before Lemma 2.5). Hence S⋆ ≤ S̄⋆ pointwise.
- Condition 1.3: S⋆(λ1,λ2) ≤ 0 for all λ1, λ2.
- Lemma 2.5 (all κ):
  - (b) S⋆(1,0) = S̄⋆(1,0) = 0;
  - (c) ∇S⋆(1,0) = ∇S̄⋆(1,0) = 0;
  - (d) ∇²S⋆(1,0) ⪯ ∇²S̄⋆(1,0).
- Inner minimisation over s:
  - ∂_s S⋆ = sψ0 - α⋆ E[N E(U+sN)].
  - Lemma 7.14: ∂²_s S⋆ > ψ0 - α⋆E N² = 0, so S⋆(Λ,·) is strictly convex.
  - Lemma 7.15: the minimiser at Λ_{1,0} is s = √(1-q0).
- **Lemma 7.17 (set O).** As printed:
  O = {Λ : d0 E[MΛ] + E[ḢΛ] > α⋆κ}.
  S⋆(Λ,s) → +∞ as s→∞ on O, and → −∞ off O, so S⋆(Λ) = −∞ off O.
  - **[note]** The proof's own display, ½α⋆E V² = ½s²ψ0 + α⋆ s E[UN] + O(1), gives a different κ term. Using α⋆E[ĤN] = q0 d0 (Gaussian integration by parts) and α⋆E N² = ψ0:
    α⋆ E[UN] = (α⋆κ **E[N]** − d0 E[MΛ] − E[ḢΛ]) / √(·).
  - So the κ-term should be α⋆κE[N]. The printed α⋆κ lacks the E[N] factor; at κ=0 the difference vanishes.
  - With the corrected form, Λ_{1,0} ∈ O for every κ: its coefficient equals −(1−q0)d0 > 0, using the identity in §4 below.

## 4. d0 and κ-general identities

- d0 = α⋆ E F'_{1-q0}(q0^{1/2}Z) < 0 (Condition 3.4).
- **[note]** By Gaussian integration by parts:
  - α⋆E[ĤN] = q0 d0 holds for all κ.
  - Eq. (109), d0 = −(1−q0)ψ0, is **κ=0 only**.
  - The general form is d0 = −(1−q0)ψ0 + α⋆κE[N]. I derived this and confirmed it numerically to 1e-14.
- Eq. (100) is κ=0 only. It reads p2(ψ0) = q0 and r2(γ0) = (1−q0)ψ0/α⋆, with r_k(γ) = E E(γ^{1/2}Z)^k.
  - Its second half uses that the κ=0 argument of E is −γ^{1/2}Z, together with symmetry.
  - For κ≠0, R uses E((κ−q^{1/2}Z)/(1−q)^{1/2}), which is not of this form.

## 5. Conditions 3.3 and 3.4 (general κ)

- **Condition 3.3.** α⋆ E[th'(ψ0^{1/2}Z)²] · E[F'_{1-q0}(q0^{1/2}Z)²] < 1, where th' = 1 − th².
  - At κ=0, Claim B.7 gives ≤ a_ub = 0.5446.
- **Condition 3.4.**
  - m(z) = E[(z + ch²(ψ0^{1/2}Z))^{-1}] on (−1,∞).
  - f̂0(x) = −F'/(1+(1−q0)F') = E'(u)/((1−q0)(1−E'(u))), with u = (κ−x)/(1−q0)^{1/2}.
  - λ(z) = z − α⋆ E[ f̂0(q0^{1/2}Z) / (1 + m(z) f̂0(q0^{1/2}Z)) ] − d0.
  - The condition is λ0 := inf_{z>−1} λ(z) < 0.
  - Lemma 3.5: λ' = 1 − α⋆θ(z), with θ(z) = E[(z+ch²)^{-2}] · E[(f̂0/(1+m f̂0))²]. θ is strictly decreasing from +∞ to 0.
  - Hence λ0 = λ(z0), where z0 = θ^{-1}(1/α⋆).
  - At κ=0, App. B uses:
    - ẑ = −0.669316
    - m(ẑ) ≤ 0.9309695
    - g(m(ẑ),q0,γ0) ≥ 0.7739, with g(m,q,γ) = E[E'(γ^{1/2}Z)/((1−q)(1−E'(γ^{1/2}Z)) + mE'(γ^{1/2}Z))] (eq. (102))
    - λ(ẑ) ≤ −0.1906
  - The κ=0 shortcut λ = ẑ − α g + (1−q0)ψ0 uses (109).

## 6. Gradient (87) and Hessian (Lemma B.9, general κ, at (1,0))

- Notation: Δ = (u1∂λ1 + u2∂λ2)Λ = (1−Λ²)(u1Ḣ + u2M).
- (87) gives ⟨∇S̄⋆,u⟩, as quoted in my notes. At (1,0) it vanishes for all κ, using the κ-correct F' identity from §1.
- **Hessian at (1,0), general κ.** This is the display in the proof of Lemma B.9 just before "Specializing further to κ=0 (which was not used up to here)".
  Set a = E[MΔ], b = E[ḢΔ], and K = (κ − Ĥ − (1−q0)N)/(1−q0). Then
  ```
  <∇²S̄⋆(1,0), u⊗u> = − E[(1−M²)(u1Ḣ+u2M)²]
      + α⋆ E{ F'_{1−q0}(Ĥ) ( −a Ĥ/q0 − b N/ψ0 + K a )² }
      − α⋆ E{ N [ (−2aĤ/q0 − 2bN/ψ0) a/(1−q0) + 3 K a²/(1−q0) + K a²/q0 ] }.
  ```
  - The paper's first term reads −E[(u1Ḣ+u2M)Δ]. At (1,0) this equals −E[(1−M²)(u·)²].
  - For u = e1: Δ = (1−M²)Ḣ. For u = e2: Δ = (1−M²)M.
  - The bilinear form follows by polarisation.
  - **[note]** My implementation of this general-κ display agrees with central finite differences of S̄⋆ to about 1e-6 for every κ tested. At κ=0 it reproduces the C1, C2, C3, I1, I2, I3 and M_ij of Claims B.10–B.12.
- **κ=0 specialisation.** With M11 = −I1 + C1I1² + C2I1I2 + C3I2², etc., and the C_i of Lemma B.9, this uses (109). It is **not valid for κ≠0**.
- Hessian of S⋆ (Lemma 2.5(d) proof):
  ∇²S⋆ = ∇²_λ S⋆(·,s0) − (∇_λ ∂_s S⋆)^{⊗2} / ∂²_s S⋆, at s0 = √(1−q0).
