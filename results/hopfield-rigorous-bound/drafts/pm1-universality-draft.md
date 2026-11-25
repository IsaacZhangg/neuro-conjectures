> **Draft, not reviewed.** Agent-written proof plan for hypothesis (U) (the ±1 step), 2026-11-25. Nothing here is claimed proved until it passes an independent review and is folded into PROOF.md. Paths to scratch scripts refer to the drafting session and are not in the repo.

# ±1 Hopfield capacity: a route that proves (U), and what it gives

Written 2026-11-25 by a research subagent. Nothing in the git repo was touched.
Companion to `/home/claude/neuro-conjectures/results/hopfield-rigorous-bound/PROOF.md`
(same notation: X is M×N with iid ±1 entries, M/N → α, Σ_k = configurations with
exactly k minus signs, k = ⌈δN⌉, δ = Q(τ), φ = normal density, ψ(t) = E(g−t)_+ =
φ(t) − tQ(t)).

## 0. Bottom line

- **Claim.** Hypothesis (U) holds, by the proof draft below, which I believe is
  complete but which nobody else has checked
  (Theorem 3 + Proposition 4 ⇒ Corollary 5 = (U)). With it, Theorem 2 of PROOF.md
  becomes unconditional: for ±1 patterns the energy barrier, and a fixed point
  within ⌈δN⌉ − 1 flips of every pattern, hold with probability → 1 for every
  α < A(1.8259) = 0.12948. That would move the rigorous ±1 bound from 0.113 to 0.1294.
- **Status.** Agent-drafted in one session, **not independently reviewed**. Every
  step is elementary or a cited standard theorem, and the one inequality that
  carries the model-specific content (Proposition 4) was also attacked
  numerically without finding a violation. Treat it as a candidate proof until a
  second reader has gone through §3. The points most likely to hide an error
  are listed in §9.
- **No bottleneck lemma is left open** for 0.1294. The step that plays that role
  is the scalar inequality (S) of Proposition 4; it is proved, with explicit
  constants.
- **Why the earlier attempts failed and this does not.** The earlier attempts ran
  Lindeberg on the whole soft maximum and needed the Gibbs measure to be
  delocalized. Here the possibly localized part is isolated deterministically:
  at any σ, at most r rows μ can have |x_μ·σ| above ‖Xσ‖/√r. Those r = εM rows
  are (i) chosen by a union bound whose entropy M·H(ε) is negligible, (ii) never
  replaced by Gaussians, and (iii) controlled by Hoeffding, not by universality.
  The gain a configuration gets from aligning with them is paid for by the
  bulk term √(NM) it gives up, and that trade is a one-dimensional inequality.
  This is the "dissection" of Oymak–Tropp (2018) applied to the left vector y,
  plus a new last step, because here the non-replaced rows are not lower order.

## 1. Statements

For ρ ∈ [0,1] and a unit vector a ∈ ℝⁿ (any n), let
Z_{ρ,a} = ρ g₀ + √(1−ρ²) Σ_j a_j ε_j, with g₀ standard normal and ε_j independent
signs. For a random variable Z put Ψ_δ(Z) = inf_t { δt + E(Z − t)_+ } (δ times the
upper δ-superquantile). Define

  V(α, δ) = sup over ρ ∈ [0,1], n ≥ 1, a ∈ S^{n−1} of { ρ√α + 2 Ψ_δ(Z_{ρ,a}) }.

**Theorem 3 (one-sided universality with one extra direction).** Let δ ∈ (0, 1/2),
k = ⌈δN⌉, M/N → α ∈ (0, ∞), X an M×N matrix of iid uniform signs. Then

  limsup_N (1/N) E max_{σ∈Σ_k} ‖Xσ‖ ≤ V(α, δ).

At ρ = 1 the bracket is √α + 2φ(τ) (the Gaussian Sudakov–Fernique value), so
V ≥ √α + 2φ(τ) always; the question is whether ρ < 1 can do better.

**Proposition 4 (scalar inequality (S)).** Let τ > 0, δ = Q(τ),
D(τ) = δ√(2 ln(1/δ)) − φ(τ) ≥ 0. If some ρ₀ ∈ (0,1) satisfies

  (a) D(τ) ≤ (1−ρ₀)√α/2  and  (b) φ(0)(1−ρ₀)(1+ρ₀)²/(6ρ₀³) ≤ √α/2,

then V(α, δ) = √α + 2φ(τ). For τ₀ = 1.8259 and ρ₀ = 0.85: D(τ₀) = 0.0129403,
(a) holds iff α ≥ 0.02977, (b) has left side 0.05558 and holds iff α ≥ 0.01236.
So (S) holds at τ₀ for every α ≥ 0.03.

**Corollary 5 ((U)).** For τ, α as in Proposition 4,
E max_{σ∈Σ_k} ‖Xσ‖ ≤ √(NM) + 2Nφ(τ) + o(N). In particular (U) holds at
τ₀ = 1.8259 for every α ∈ [0.03, A(τ₀)), A(τ₀) = 0.129489….

**Consequence (Theorem 2 of PROOF.md without the hypothesis).** For ±1 patterns
and α ∈ [0.03, 0.12948): with probability 1 − O(N e^{−cN}), every stored pattern
has an energy barrier at radius ⌈δ(τ₀)N⌉ and a fixed point (zero-field
convention) inside it, overlap ≥ 0.932. For α < 0.03 use another τ (a
floating-point grid search finds, for every α ∈ [10⁻⁶, 0.12948], a τ with
α < A(τ) and (a),(b); not certified) or the classical bounds (Newman 0.056).

## 2. Proof of Proposition 4 (proved; elementary)

Write c = √(1−ρ²), s = Σ a_j ε_j, Z = ρg₀ + cs.

1. *Value at ρ = 1.* d/dt [δt + ψ(t)] = δ − Q(t) vanishes at t = τ, so
   Ψ_δ(g₀) = δτ + φ(τ) − τQ(τ) = φ(τ).
2. *Hoeffding bound, any ρ, a.* E e^{λZ} = e^{λ²ρ²/2} Π_j cosh(λca_j) ≤ e^{λ²/2}.
   Since x_+ ≤ e^{λx}/(eλ), E(Z−t)_+ ≤ e^{λ²/2−λt}/(eλ). Choosing t with
   e^{λ²/2−λt−1} = δ gives δt + E(Z−t)_+ ≤ δ(λ/2 + ln(1/δ)/λ); λ = √(2 ln(1/δ))
   gives Ψ_δ(Z) ≤ δ√(2 ln(1/δ)) = φ(τ) + D(τ).
3. *Fourth-moment bound, ρ > 0.* Let F(x) = E(ρg₀ + x − τ)_+ = ρψ((τ−x)/ρ). Then
   E(Z−τ)_+ = E F(cs), and E F(cg′) = ψ(τ) for g′ standard normal. F is smooth
   with F''''(x) = ρ⁻³ φ''((τ−x)/ρ), and sup_u |φ''(u)| = sup |u²−1|φ(u) = φ(0).
   Replace a_jε_j by a_jg_j one j at a time. In the Taylor expansion of
   F(x + bξ) to fourth order, orders 0–3 agree for ξ = ε and ξ = g (symmetric,
   variance 1), and the remainders are at most b⁴ sup|F''''| (Eε⁴ + Eg⁴)/24 =
   b⁴ sup|F''''|/6. With b = ca_j and Σ a_j⁴ ≤ 1:
   E(Z−τ)_+ ≤ ψ(τ) + φ(0)c⁴/(6ρ³), hence
   Ψ_δ(Z) ≤ δτ + E(Z−τ)_+ ≤ φ(τ) + φ(0)(1−ρ²)²/(6ρ³).
4. *Combine.* If ρ ≤ ρ₀: ρ√α + 2Ψ_δ ≤ √α + 2φ(τ) + 2D − (1−ρ₀)√α ≤ √α + 2φ(τ) by (a).
   If ρ ≥ ρ₀: 2φ(0)(1−ρ²)²/(6ρ³) = (1−ρ)·φ(0)(1−ρ)(1+ρ)²/(3ρ³), and
   q(ρ) = (1−ρ)(1+ρ)²/ρ³ is decreasing on (0,1] (d/dρ log q =
   −1/(1−ρ) + 2/(1+ρ) − 3/ρ < 0), so the factor is at most its value at ρ₀,
   which is ≤ √α by (b); hence ρ√α + 2Ψ_δ ≤ √α + 2φ(τ) − (1−ρ)√α + (1−ρ)√α. ∎

Remarks. The slack is large: Hoeffding covers ρ ≤ 0.928 and the fourth-moment
bound covers ρ ≥ 0.69 at α = 0.1294. The inequality is *not* universal in
(α, δ): for the whole cube (δ = 1/2) and small α it fails, and it should,
because there max ‖Xσ‖² ≥ N² + (M−1)N (take σ = a row) exceeds the Gaussian
bound N²(√α + √(2/π))² when α < 0.0518. Theorem 3 reproduces that: with
a = e₁, ρ = √α it gives 1 + α/2 + O(α²). So (U) is true for the Hamming sphere
at these parameters for a reason (thin tail δ = 0.034 plus the √α penalty), not
by a general principle.

## 3. Proof of Theorem 3 (proved; standard tools)

Fix ε ∈ (0,1), η ∈ (0,1/2), β > 0, r = ⌈εM⌉. C denotes universal constants.

**Step 1 (dissection; deterministic).** For T ⊂ [M] let
B_T = { y ∈ ℝ^M : ‖y‖ ≤ 1, |y_μ| ≤ r^{−1/2} for μ ∉ T } and
F_T(A) = max_{σ∈Σ_k, y∈B_T} yᵀAσ.

*Lemma 1.* For every matrix A, max_{σ∈Σ_k} ‖Aσ‖ = max_{|T|=r} F_T(A).

*Proof.* "≥" because B_T is inside the unit ball. "≤": take σ, u = Aσ ≠ 0, and
T = the indices of the r largest |u_μ|. The r+1 largest squares sum to at most
‖u‖², so |u_μ| ≤ ‖u‖/√(r+1) for μ ∉ T. Thus y = u/‖u‖ ∈ B_T and yᵀAσ = ‖u‖. ∎

**Step 2 (union over T).** F_T is a maximum of linear functions of the entries
with coefficient vectors yσᵀ of Frobenius norm ≤ √N: convex and √N-Lipschitz.
Talagrand's inequality (Ledoux, Cor. 4.10, entries in [−1,1]) gives
P(|F_T(X) − med| ≥ s) ≤ 4e^{−s²/(16N)}. All F_T with |T| = r have the same law
(rows are exchangeable). Integrating the union bound over J = binom(M, r) sets,

  E max_σ ‖Xσ‖ ≤ E F_{T₀}(X) + 4√(N ln(4J)) + 17√N,  T₀ = {1,…,r},

and ln J ≤ r ln(eM/r) ≤ (εM + 1) ln(e/ε).

**Step 3 (Lindeberg on the capped rows only).** Let H be X with the rows outside
T₀ replaced by independent standard Gaussian rows. Let Y ⊂ B_{T₀} be an η-net
of B_{T₀}, |Y| ≤ (1+2/η)^M, K = Σ_k × Y, and
f_β(A) = β⁻¹ ln Σ_{(σ,y)∈K} exp(β yᵀAσ). Then, for any A,

- max_K ≤ f_β(A) ≤ max_K + β⁻¹(N ln 2 + M ln(1+2/η)),
- max_K ≤ F_{T₀}(A) ≤ max_K + η‖A‖_op √N,
- for μ ∉ T₀: ∂³f_β/∂A_μi³ = β² κ₃(y_μσ_i) (third cumulant under the Gibbs
  measure), and |y_μσ_i| ≤ r^{−1/2}, so |∂³f_β/∂A_μi³| ≤ 2β² r^{−3/2}, uniformly
  in all other entries.

Swap the N(M−r) entries outside T₀ one at a time (first two moments match):
each swap costs at most (2β²r^{−3/2}/6)(E|ε|³ + E|g|³) ≤ β² r^{−3/2}. With
E‖X‖_op ≤ √(π/2)(√N + √M) (contraction plus Gordon/Chevet),

  E F_{T₀}(X) ≤ E F_{T₀}(H) + Cη(1+√α)N + β⁻¹(ln 2 + (M/N) ln(1+2/η))N + β² N/(ε^{3/2}√M).

For fixed ε, η, β the last term is O(√N). (This step is also Han–Shen 2023,
Theorem 2.5, applied conditionally on X_{T₀} with u = √M·y_{T₀ᶜ},
L_u = ε^{−1/2}, L_w = 1 and Q(u,w) = y_{T₀}ᵀX_{T₀}σ.)

**Step 4 (Sudakov–Fernique, conditionally on X_{T₀}).** Write T = T₀, G for the
Gaussian rows, y = (y_T, y′). Given X_T, compare on Σ_k × B_T

  𝒳(σ,y) = y′ᵀGσ + y_TᵀX_Tσ  with  𝒴(σ,y) = √N g·y′ + ‖y′‖ h·σ + y_TᵀX_Tσ,

g ∈ ℝ^{M−r}, h ∈ ℝ^N standard Gaussian. The means agree, and with ρ = ‖y′‖,
ρ̃ = ‖ỹ′‖,

  E(𝒴−𝒴̃)² − E(𝒳−𝒳̃)² = N(ρ−ρ̃)² + 2(ρρ̃ − y′·ỹ′)(N − σ·σ̃) ≥ 0

(expand: E(𝒳−𝒳̃)² = N(ρ²+ρ̃²) − 2(y′·ỹ′)(σ·σ̃);
E(𝒴−𝒴̃)² = 2N(ρ²+ρ̃²) − 2N y′·ỹ′ − 2ρρ̃ σ·σ̃). Sudakov–Fernique with equal
means (Vitale 2000; Chatterjee 2005), via finite subsets as in Lemma 4 of
PROOF.md, then dropping the caps and maximizing over y:

  E F_T(H) ≤ E max_{σ∈Σ_k} max_{θ∈[0,π/2]} [ cos θ (√N‖g‖ + h·σ) + sin θ ‖X_Tσ‖ ].

**Step 5 (nets on the r uncapped rows).** ‖X_Tσ‖ ≤ (1−η)⁻¹ max_{a∈𝒜} aᵀX_Tσ for an
η-net 𝒜 of S^{r−1}, |𝒜| ≤ (1+2/η)^r; take an η-net Θ of [0,π/2]. For
(a,θ) ∈ 𝒜×Θ put w = X_Tᵀa ∈ ℝ^N and

  Φ_{a,θ} = cos θ √N‖g‖ + max_{σ∈Σ_k} (cos θ h + sin θ w)·σ.

The discretization error is at most Cη(√N‖g‖ + √N‖h‖ + ‖X_T‖_op√N), of mean
≤ Cη(1+√α)N. Each Φ_{a,θ} is √(2N)-Lipschitz in (g,h) and its conditional mean
is convex and √N-Lipschitz in X_T, so by Gaussian concentration and Talagrand
it has sub-Gaussian fluctuations of order √N; over |𝒜||Θ| points,

  E F_T(H) ≤ max_{a,θ} E Φ_{a,θ} + Cη(1+√α)N + C√(N(r ln(1+2/η) + ln(2/η) + 1)).

**Step 6 (one fixed direction is an iid problem).** For fixed unit a, the
coordinates w_i = Σ_{μ∈T} a_μ X_μi are iid copies of Σ a_μ ε_μ (columns of X_T
are iid). So Z_i = cos θ h_i + sin θ w_i are iid copies of Z_{cos θ, a}, symmetric.
For σ ∈ Σ_k with minus set S and any t ≥ 0,
Z·σ = Σ_i Z_i + 2Σ_{i∈S}(−Z_i) ≤ Σ_i Z_i + 2[kt + Σ_i(−Z_i − t)_+], hence, with
E‖g‖ ≤ √M, k ≤ δN + 1, and the minimizing t ≤ √(2 ln(1/δ)),

  E Φ_{a,θ} ≤ N [ cos θ √(M/N) + 2Ψ_δ(Z_{cos θ,a}) ] + 2√(2 ln(1/δ)) ≤ N·V(M/N, δ) + C_δ.

**Conclusion.** Chaining Steps 1–6 and letting N → ∞:

  limsup (1/N) E max ‖Xσ‖ ≤ V(α,δ) + C√(αε ln(e/ε)) + Cη(1+√α)
    + (ln 2 + α ln(1+2/η))/β + C√(αε ln(1+2/η)).

Let β → ∞, then ε → 0, then η → 0. (V is 1-Lipschitz in √α.) ∎

Remarks.
- Steps 1–5 use nothing about Σ_k: for any deterministic Σ ⊂ {±1}^N,
  E max_Σ ‖Xσ‖ ≤ sup_{θ,a} E[cos θ √(NM) + max_Σ (cos θ h + sin θ w_a)·σ] + o(N).
  This is the zero-temperature analogue of "Gaussian patterns plus finitely many
  Bernoulli fields" (Ali et al. 2026; Barra–Durante 2026).
- The proof is qualitative. The rates are terrible (ε must be tiny before the
  union bound is small, and the Lindeberg term has ε^{−3/2}); this does not
  matter for a limit statement.
- Only Hoeffding and matching first two moments are used, so the same proof
  covers any symmetric entries with variance 1 and variance proxy 1.

## 4. Evaluation of the proposed ideas

1. **Row-wise one-sided replacement at zero temperature.** As a stand-alone
   comparison E f(Rademacher rows) ≤ E f(Gaussian rows): **false in general**
   (whole cube, α < 0.0518; and signs are not below Gaussians in convex
   order). Its kernel is right and is used: Hoeffding for the few rows that
   are not replaced (Proposition 4, item 2). Hoeffding alone on those rows
   gives only V ≤ √α + 2δ√(2 ln(1/δ)), i.e. α ≈ 0.0880. The fourth-moment
   comparison near ρ = 1 (item 3) is what recovers the full constant.
2. **Han–Shen / Montanari–Saeed.** Checked in the papers. Han–Shen Theorem 2.5
   needs both index sets inside ℓ∞ boxes, S_u ⊂ [−L_u,L_u]^m, S_w ⊂ [−L_w,L_w]^n,
   with error polynomial in L; it applies to the capped problem (L_u = ε^{−1/2})
   and can replace Step 3 as a black box. Montanari–Saeed Theorem 1 needs a
   nonnegative loss that is Lipschitz or locally Lipschitz with two-sided
   polynomial growth, and a parameter set inside a domain of Gaussianity; our
   loss −u² is neither, and it cannot be covered, since universality of
   max Σ(x_μ·σ)² over a delocalized set is false (same counterexample).
   **Verdict: the capped-y reduction is the right first half; the
   large-coordinate part is not negligible and needs Steps 2, 4–6.** This is
   the route taken.
3. **Cheaper bounds above 0.113.** None found. Splitting into bulk + worst-case
   tail gives 0.1013 (the worst-case tail is about 4φ²N², not 4δ²N²: εM rows can be
   aligned in a Gaussian-like direction a, giving ‖X_Tσ‖ ≈ 2Nφ). Hoeffding-only gives 0.0880.
   Contraction gives 0.0282. The joint treatment of bulk and tail is essential.
4. **Reference state (0.1347).** See §6: it transfers to ±1, conditionally on the
   Gaussian-side statement.

## 5. Constants (floating point unless marked)

| Route | ±1 bound | Status |
| :-- | :-- | :-- |
| Theorem 3 + Proposition 4 ⇒ (U) ⇒ Theorem 2 | 0.12948 at τ = 1.8259 (constant certified by check.py) | proved here, unreviewed |
| Same plus reference state σ₁ (§6) | 0.13472 at τ = 1.942, θ′ = 2.531 | plausible; needs the Gaussian-side lemma |
| Bulk universality + worst-case tail (sum of maxima) | 0.1013 at τ = 2.37 | would be provable; below 0.113 |
| Hoeffding only on uncapped rows | 0.0880 at τ = 2.22 | provable; below 0.113 |
| Contraction, factor √(π/2) | 0.0282 | provable; useless |
| Annealed union bound | 0.056 (Newman) | known |

## 6. Reference-state gain: does it transfer to ±1?

Reconstruction of the Gaussian heuristic (not proved here). Flip the spins
with noise field below −θ, θ′ = θ/√α, p = Q(θ′). Then
‖Gσ₁‖²/N² → (√α + 2φ(θ′))² − 8Q(θ′)√α φ(θ′) =: g_ref(θ′), and the criterion is

  (√α + 2φ(τ))² + (1−2Q(τ))² < g_ref(θ′) + (1−2Q(θ′))²  for some θ′ > τ.

This gives α* = 0.13472, matching the reported ≈ 0.1347.

Transfer to ±1 needs two things.

- *Upper side.* (U) at the new radius: Proposition 4 at τ = 1.942 has
  D = 0.009875 and holds (ρ₀ = 0.85) for all α ≥ 0.0174. **Covered.**
- *Lower side.* Use as reference the best state at radius k′ = ⌈pN⌉ (Lemma 2
  works with any reference inside the ball). max_{Σ_{k′}} ‖Xσ‖ is convex
  Lipschitz, so it concentrates, and a lower bound on its mean suffices. With
  the capped ball B^L = {‖y‖ ≤ 1, ‖y‖_∞ ≤ L/√M}, Step 3 with T = ∅ is two-sided:

  **Proposition 6 (lower-bound transfer; proved by Step 3).**
  E max_{σ∈Σ} ‖Xσ‖ ≥ E max_{σ∈Σ, y∈B^L} yᵀGσ − o(N) for every fixed L and every
  deterministic Σ ⊂ {±1}^N.

  So any Gaussian lower bound whose witness y is delocalized transfers. For
  σ₁ the witness y ∝ Gσ₁ has Gaussian-like coordinates, and truncating at
  L√N loses o_L(1). **What is missing is on the Gaussian side:** a proof that
  E max_{σ∈Σ_{k′}, y∈B^L} yᵀGσ ≥ N(√g_ref(θ′) − o_L(1)), i.e. the one-step
  computation made rigorous, with exactly k′ flips and a capped witness.

## 7. Numerical checks (scripts in `/tmp/claude-0/-home-user/7b4273ff-69b1-55a5-830b-99eb38aec9c1/scratchpad/work/`)

- `scalar_check.py`: exact evaluation of ρ√α + 2Ψ_δ(Z_{ρ,a}) − √α − 2φ(τ) on a ρ-grid
  for equal weights n = 1…40, 60, 100, 200 and 300 random weight vectors
  (n ≤ 12). Maximum is 0, attained only at ρ = 1. Also sup_{n<400} Ψ_δ(S_n/√n) =
  0.075326 < φ(τ) = 0.075333 (approached from below).
- `scalar_adversarial.py`: Nelder–Mead over (ρ, a), n ≤ 8, 40 restarts, for
  (α, τ) = (0.1294, 1.8259), (0.05, ·), (0.02, ·), (0.005, ·), (0.1347, 1.942).
  No positive excess found.
- `prop4_numbers.py`, `grid_cert.py`: the constants in Proposition 4 (30 digits),
  and the grid search over α ∈ [10⁻⁶, 0.12948].
- `mc_check.py`: exhaustive maxima at N = 24–40, 400 trials. The exact hybrid
  inequality of Step 4 holds within noise (e.g. N = 40, M = 5, k = 3, r = 1:
  23.93 ± 0.20 ≤ 24.26 ± 0.33). Rademacher and Gaussian maxima are within
  noise of each other and below the SF bound. Finite-N evidence only.
- `constants.py`: the table in §5 and α* = 0.13472.

None of this is a proof. The certified part of the final constant is still
check.py's A(1.8259).

## 8. Literature checked (through the paper index; arxiv.org is blocked)

- Han, Shen, arXiv:2206.07936, Theorem 2.5: read the statement; ℓ∞ boxes on
  both sides.
- Montanari, Saeed, arXiv:2202.08832, Assumptions 2–5 and Theorem 1: read.
- Oymak, Tropp, arXiv:1511.09433, §9.3–9.4: dissection T_J (coordinates outside
  J at most (#J)^{−1/2}), binom(n,k) pieces, partial Lindeberg replacement.
  Same device as Steps 1–3, used there for a *lower* bound on a minimum, where
  the non-replaced part is lower order.
- Stojnic, arXiv:2403.01907: the binary switch is one sentence ("concentrations
  and the Lindeberg variant of the central limit theorem allow for such a
  switch"), no proof.
- Barra–Durante arXiv:2610.41938 and Ali et al. arXiv:2611.03306 (abstracts):
  positive-temperature results with the same structure (Gaussian bulk plus
  finitely many binary patterns or Bernoulli fields).
- No paper found that proves (U) or a ±1 bound above 0.113.
- From memory, not re-read: Vitale 2000 and Chatterjee 2005 (Sudakov–Fernique
  with equal means), Ledoux Cor. 4.10, Chatterjee 2006 (Lindeberg).

## 9. What a reviewer should attack; go / no-go

1. Step 4: the increment identity, and that Sudakov–Fernique is used with a
   common non-zero mean (if in doubt, rerun the interpolation proof; only
   covariances enter).
2. Step 3: the third-derivative bound uses the cap only for the rows that are
   swapped; rows in T₀ are never differentiated.
3. Steps 2 and 5: the two union bounds (binom(M,r) sets; (1+2/η)^r directions)
   and the order of limits β → ∞, ε → 0, η → 0.
4. Proposition 4, item 3: the fourth-order remainder and sup|φ''| = φ(0).
5. Step 6: that w has iid coordinates for a fixed a (it does: columns of X_T).
6. Whether "all α < 0.12948" should be stated via the grid search (uncertified)
   or as α ∈ [0.03, 0.12948) plus classical results below.

**Go.** Recommended next actions: independent line-by-line review of §2–3; then
fold Theorem 3, Proposition 4 and Corollary 5 into PROOF.md as §6, restate
Theorem 2 without (U), and extend check.py with the three numbers of
Proposition 4 (D(τ₀), and the two thresholds) in exact arithmetic. The 0.1347
extension is a separate, Gaussian-side task (§6).
