> **Not yet independently reviewed.** Agent-written proof of the ±1 universality step, 2026-11-25. The result README stays Partial until an independent review passes.

# ±1 Hopfield model: an unconditional proof of Theorem 2

Companion to `/home/claude/neuro-conjectures/results/hopfield-rigorous-bound/PROOF.md`
(cited below as PROOF.md). Same notation. This file proves Theorem 2 of PROOF.md
without hypothesis (U), for every load α ∈ (0, 0.12948), with one radius for all α.

## 0. Status

- **What is proved.** Theorem 2′ below: for iid uniform ±1 patterns with P/N → α ∈
  (0, 0.12948), with probability 1 − O(N e^{−cN}) *every* stored pattern ξ^ν has an
  energy barrier on the Hamming sphere of radius k = ⌈δ₀N⌉, δ₀ = Q(1.8259) ≈ 0.03393,
  and hence a fixed point σ* (σ*_i h_i(σ*) ≥ 0 for all i, the zero-field-keeps
  convention of PROOF.md Lemma 2′; strict if there are no ties) with
  d(σ*, ξ^ν) ≤ k − 1, overlap ≥ 0.932.
- **What it rests on.** Every step is written out here except four standard results
  that are quoted: Talagrand's concentration inequality for convex Lipschitz
  functions of bounded independent variables, Gaussian concentration for Lipschitz
  functions, Gaussian integration by parts, and the volume bound for nets. Their
  theorem numbers were taken from memory and not re-checked against the books.
  All numbers are certified in exact rational arithmetic by `check_pm1.py` (§8).
- **Nothing is left conditional.** No step had to be labelled conditional.
- **Review status.** Written by an agent from a planner's draft, checked
  adversarially in one session (§9 lists what changed), not independently reviewed
  by a person. Treat it as a candidate proof until a second reader has gone through
  §3–§7.

## 1. Setting and statement

Patterns ξ^1, …, ξ^P ∈ {±1}^N are iid uniform. J_ij = (1/N) Σ_ν ξ^ν_i ξ^ν_j for i ≠ j,
J_ii = 0, h_i(σ) = Σ_{j≠i} J_ij σ_j, H(σ) = −½ Σ_{i≠j} J_ij σ_i σ_j, d = Hamming distance.

Constants (all certified in §8):

| symbol | definition | value |
| :-- | :-- | :-- |
| τ₀ | 18259/10000 | 1.8259 |
| δ₀ | Q(τ₀), Q = standard normal upper tail | 0.0339326631… |
| φ₀ | φ(τ₀), φ = standard normal density | 0.0753287102… |
| L₀ | ln(1/δ₀) | 3.3833772150… |
| φ_H | δ₀ √(2L₀) | 0.0882689891… |
| D | φ_H − φ₀ | 0.0129402789… |
| t_H | (2L₀ − 1)/√(2L₀) | ∈ [0, √(2L₀)] |
| t* | √(2L₀) | ≥ τ₀ (C1) |
| ρ₀ | 17/20 | 0.85 |
| A(τ₀) | ((δ₀(1−δ₀) − φ₀²)/φ₀)² | 0.1294899826… |
| A_H | (1 − δ₀ − 2δ₀L₀)²/(2L₀) | 0.0801512047… |

k = k(N) = ⌈δ₀N⌉, Σ_k = {σ ∈ {±1}^N with exactly k entries −1}, **1** = all-ones vector.

**Theorem 2′.** Let P = P(N) with P/N → α ∈ (0, 0.12948). There are c > 0 and N₀
such that for N ≥ N₀, with probability at least 1 − 5P e^{−cN}, for every
ν ∈ {1, …, P}:

1. H(σ) > H(ξ^ν) for every σ with d(σ, ξ^ν) = k;
2. there is σ* with d(σ*, ξ^ν) ≤ k − 1 and σ*_i h_i(σ*) ≥ 0 for every i (and
   σ*_i h_i(σ*) > 0 for every i whenever H takes different values at any two
   configurations that differ in one spin, as in PROOF.md Lemma 2); in particular
   ξ^ν·σ*/N ≥ 1 − 2δ₀ ≥ 0.932.

Since P ≤ 2αN for large N, 5P e^{−cN} → 0. The proof splits into Case A,
α ∈ (3/100, 0.12948), and Case B, α ∈ (0, 3/100]. Both use the same δ₀.

## 2. Deterministic reductions

**Gauge.** Fix ν and write ξ = ξ^ν. Put M = P − 1 and let X = X^{(ν)} be the M × N
matrix with rows x_μ = (ξ^μ_i ξ_i)_i over μ ≠ ν. Given ξ, the entries ξ^μ_i ξ_i are
iid uniform signs, so X has iid uniform ±1 entries (and is independent of ξ). The map
σ ↦ (ξ_i σ_i)_i is a bijection of {±1}^N that preserves Hamming distances, maps ξ to
**1**, and (since σ_i h_i(σ) and H are invariant under simultaneously gauging σ and
the patterns) it is enough to prove 1–2 for the gauged system, where by PROOF.md
Lemma 1 (with Σ_{μ,i} X_μi² = MN a constant)

  H(σ) = −(1/(2N)) [ (**1**·σ)² + ‖Xσ‖² ] + const.

**Barrier criterion.** For σ ∈ Σ_k, **1**·σ = N − 2k, so

  H(σ) > H(**1**)  ⟺  ‖Xσ‖² < ‖X**1**‖² + 4k(N − k).

Let F := max_{σ∈Σ_k} ‖Xσ‖. Statement 1 holds as soon as

  (★)  F² < ‖X**1**‖² + 4k(N − k).

**Fixed point.** By PROOF.md Lemma 2′ (proof: let σ* minimize H on the ball
B = {d(σ, **1**) ≤ k − 1}; a single flip either stays in B, so cannot lower H, or lands
on the sphere, where H > H(**1**) ≥ H(σ*); the flip changes H by 2σ*_i h_i(σ*)),
statement 1 implies statement 2 with weak inequalities, and Lemma 2 of PROOF.md
gives strict ones when there are no ties. Ties can occur for ±1 patterns (energies
are rationals with denominator 2N), which is why the weak form is the statement.
Finally d(σ*, **1**) ≤ k − 1 ≤ δ₀N gives overlap ≥ 1 − 2δ₀ ≥ 0.932 (C8).

So Theorem 2′ follows from: P((★) fails) ≤ 5e^{−cN} for each fixed ν (the law of
X^{(ν)} does not depend on ν), and a union bound over ν. The work is an upper bound
on F (§3–§6) and a lower bound on ‖X**1**‖² (§7).

## 3. Tools

**T1 (Talagrand).** Let W ∈ ℝ^n have independent coordinates with values in
[−1, 1], and let f: ℝ^n → ℝ be convex and L-Lipschitz for the Euclidean norm. Then
for every median m_f of f(W) and every s ≥ 0,

  P(|f(W) − m_f| ≥ s) ≤ 4 exp(−s²/(16L²)).

*Source.* Talagrand (1988); Ledoux (2001), Corollary 4.10, states
P(|f − m| ≥ s) ≤ 4e^{−s²/4} for convex 1-Lipschitz f on [0, 1]^n. With W = 2U − 1,
U ∈ [0, 1]^n, the function U ↦ f(2U − 1)/(2L) is convex and 1-Lipschitz, which gives
the constant 16.

**Lemma 3.1 (median versus mean).** Under T1, |E f(W) − m_f| ≤ 8L√π.

*Proof.* |E f − m_f| ≤ E|f − m_f| = ∫₀^∞ P(|f − m_f| ≥ s) ds ≤ 4 ∫₀^∞ e^{−s²/(16L²)} ds
= 4 · (√(16πL²))/2 = 8L√π. ∎

With L = √N this is 8√(πN). Consequently, under T1 with L = √N,

  P(f ≥ E f + 8√(πN) + s) ≤ P(f ≥ m_f + s) ≤ 4e^{−s²/(16N)}.          (3.1)

**T2 (Gaussian concentration).** If W is a standard Gaussian vector in ℝ^n and f is
L-Lipschitz, then P(f(W) − E f(W) ≥ s) ≤ e^{−s²/(2L²)} (Tsirelson–Ibragimov–Sudakov;
Boucheron, Lugosi, Massart (2013), Theorem 5.6).

**Lemma 3.2 (union-bound integration).** Let W_1, …, W_J be random variables with
P(W_j ≥ s) ≤ C e^{−s²/v} for all s ≥ 0 and j, where C ≥ 1, v > 0, CJ > 1. Then

  E max_j (W_j)_+ ≤ √(v ln(CJ)) + √v / (2√(ln(CJ))).

*Proof.* Let u = √(v ln(CJ)). For any u ≥ 0,
E max_j (W_j)_+ ≤ u + ∫_u^∞ P(max_j W_j ≥ s) ds ≤ u + CJ ∫_u^∞ e^{−s²/v} ds. By the
Mills bound ∫_u^∞ e^{−s²/v} ds ≤ (v/(2u)) e^{−u²/v} = v/(2uCJ). So the right side is
at most u + v/(2u). ∎

With v = 16N the second term is 2√N/√(ln(CJ)) ≤ 2√N/√(ln C).

**Lemma 3.3 (nets).** Let U be a subset of the closed unit ball of ℝ^n and
η ∈ (0, 1]. There is Y ⊂ U with |Y| ≤ (1 + 2/η)^n such that every u ∈ U is within
Euclidean distance η of Y.

*Proof.* Take Y ⊂ U maximal with pairwise distances > η. Maximality gives the
covering property. The balls of radius η/2 around points of Y are disjoint and lie in
the ball of radius 1 + η/2, so |Y| (η/2)^n ≤ (1 + η/2)^n. ∎

**Lemma 3.4 (operator norm).** Let W be an m × n matrix with independent centered
entries satisfying E e^{λW_ab} ≤ e^{λ²/2} for all real λ (true for uniform signs,
since cosh x ≤ e^{x²/2}, and for standard Gaussians). Then
E‖W‖_op ≤ (21/5)√(m + n) + 2.

*Proof.* Let 𝒩_m, 𝒩_n be (1/4)-nets of the unit spheres S^{m−1}, S^{n−1} from Lemma 3.3,
of sizes ≤ 9^m, 9^n. If u*, v* attain ‖W‖_op = u*ᵀWv* and u₀, v₀ are net points within
1/4, then u*ᵀWv* = u₀ᵀWv₀ + (u* − u₀)ᵀWv* + u₀ᵀW(v* − v₀) ≤ Z + ‖W‖_op/2 with
Z = max over the nets of uᵀWv; so ‖W‖_op ≤ 2Z. For fixed unit u, v,
E e^{λuᵀWv} = Π E e^{λu_a v_b W_ab} ≤ e^{λ²/2}, hence P(uᵀWv ≥ t) ≤ e^{−t²/2} and
P(Z ≥ t) ≤ 9^{m+n} e^{−t²/2}. With t₀ = √(2(m+n) ln 9) ≥ √(4 ln 9) > 1,
E Z_+ ≤ t₀ + 9^{m+n} ∫_{t₀}^∞ e^{−t²/2} dt ≤ t₀ + 1/t₀ ≤ t₀ + 1. So
E‖W‖_op ≤ 2√(2 ln 9) √(m+n) + 2 ≤ (21/5)√(m+n) + 2 (C8). ∎

**Theorem 3.5 (Sudakov–Fernique with a common mean).** Let S be finite, let
(X_s)_{s∈S} and (Y_s)_{s∈S} be centered Gaussian vectors, and m ∈ ℝ^S deterministic.
If E(X_s − X_t)² ≤ E(Y_s − Y_t)² for all s, t ∈ S, then

  E max_s (m_s + X_s) ≤ E max_s (m_s + Y_s).

*Proof.* Only the two laws matter, so take X and Y independent, with covariances
C^X, C^Y, and put Δ = C^Y − C^X, γ^W_{st} = E(W_s − W_t)² = C^W_ss + C^W_tt − 2C^W_st.
For β > 0 let F_β(z) = β^{−1} ln Σ_s e^{βz_s}; then max_s z_s ≤ F_β(z) ≤ max_s z_s + ln|S|/β,
∂_s F_β = p_s := e^{βz_s}/Σ_t e^{βz_t}, and ∂_s∂_t F_β = β(p_s 1_{s=t} − p_s p_t); all
derivatives are bounded. Put Z(u) = √u Y + √(1−u) X + m and I(u) = E F_β(Z(u)),
u ∈ [0, 1]. I is continuous on [0, 1] by dominated convergence (|F_β(z)| ≤
max_s |z_s| + ln|S|/β). For u ∈ (0, 1), differentiating under E (the derivative of the
integrand is bounded by Σ_s (|Y_s|/(2√u) + |X_s|/(2√(1−u))), integrable uniformly on
compact subsets of (0, 1)),

  I′(u) = Σ_s E[ p_s(Z(u)) ( Y_s/(2√u) − X_s/(2√(1−u)) ) ].

Gaussian integration by parts (E W_s f(W) = Σ_t C^W_st E ∂_t f(W) for a centered
Gaussian vector W and C¹ f with bounded derivatives), applied conditionally on X to
f(y) = p_s(√u y + √(1−u)X + m) and conditionally on Y to the X-term, gives
E[Y_s p_s(Z)] = √u Σ_t C^Y_st E ∂_s∂_t F_β(Z) and E[X_s p_s(Z)] = √(1−u) Σ_t C^X_st E ∂_s∂_t F_β(Z). Hence

  I′(u) = ½ Σ_{s,t} Δ_st E ∂_s∂_t F_β(Z(u)) = (β/2) E[ Σ_s p_s Δ_ss − Σ_{s,t} p_s p_t Δ_st ].

Because Σ_t p_t = 1, Σ_s p_s Δ_ss = ½ Σ_{s,t} p_s p_t (Δ_ss + Δ_tt), so the bracket equals
½ Σ_{s,t} p_s p_t (Δ_ss + Δ_tt − 2Δ_st) = ½ Σ_{s,t} p_s p_t (γ^Y_st − γ^X_st) ≥ 0. So I is
non-decreasing and E max(m + X) ≤ I(0) ≤ I(1) ≤ E max(m + Y) + ln|S|/β. Let β → ∞. ∎

The mean m enters only as the point at which F_β is evaluated, which is why a common
non-zero mean is allowed. (This is the interpolation proof of Chatterjee (2005);
Vitale (2000) also treats equal means.)

**Corollary 3.6 (compact index sets).** Let S = S₀ × U with S₀ finite and U ⊂ ℝ^d
compact, and let X, Y, m be as in Theorem 3.5 on every finite subset of S, with
s ↦ m_s + X_s and s ↦ m_s + Y_s a.s. continuous in the U-coordinate, and
E sup_S |m + X| < ∞. Then E sup_S (m + X) ≤ E sup_S (m + Y).

*Proof.* Take a countable dense D ⊂ U and finite D_n ↑ D. By continuity
max_{S₀×D_n}(m + X) ↑ sup_S (m + X), and the first term is integrable, so by monotone
convergence and Theorem 3.5, E sup_S (m + X) = lim_n E max_{S₀×D_n}(m + X) ≤
lim sup_n E max_{S₀×D_n}(m + Y) ≤ E sup_S (m + Y). ∎

## 4. The scalar inequalities

For r ≥ 1, a unit vector a ∈ S^{r−1} and ρ ∈ [0, 1] let

  Z_{ρ,a} = ρ g₀ + √(1−ρ²) Σ_{j=1}^r a_j ε_j,

with g₀ standard normal and ε_j independent uniform signs, all independent. Z_{ρ,a} is
symmetric. Let ψ(u) = E(g₀ − u)_+ = φ(u) − uQ(u) (so ψ′ = −Q, ψ″ = φ) and

  Γ(ρ, a, t; α′) = ρ√α′ + 2( δ₀ t + E(Z_{ρ,a} − t)_+ ).

**Lemma 4.1 (value at the Gaussian point).** δ₀τ₀ + ψ(τ₀) = φ₀.

*Proof.* δ₀τ₀ + φ(τ₀) − τ₀Q(τ₀) = φ₀ because δ₀ = Q(τ₀). ∎

**Lemma 4.2 (Hoeffding).** For all r, a, ρ: δ₀ t_H + E(Z_{ρ,a} − t_H)_+ ≤ φ_H.

*Proof.* Write c = √(1−ρ²). For λ > 0, E e^{λZ} = e^{λ²ρ²/2} Π_j cosh(λ c a_j) ≤
e^{λ²ρ²/2 + λ²c²Σa_j²/2} = e^{λ²/2}. Since x_+ ≤ e^{λx}/(eλ) for all real x (the
maximum of x e^{−λx} over x > 0 is 1/(eλ)), E(Z − t)_+ ≤ e^{λ²/2 − λt − 1}/λ. Take
λ = √(2L₀) and t = t_H = λ/2 + (L₀ − 1)/λ = (2L₀ − 1)/√(2L₀); then
e^{λ²/2 − λt − 1} = e^{−L₀} = δ₀ and δ₀t + E(Z − t)_+ ≤ δ₀(t + 1/λ) = δ₀(λ/2 + L₀/λ) =
δ₀√(2L₀) = φ_H. Also 0 ≤ t_H ≤ √(2L₀) because L₀ > 1/2 (C1). ∎

**Lemma 4.3 (fourth-moment comparison).** For ρ ∈ (0, 1] and all r, a:

  δ₀τ₀ + E(Z_{ρ,a} − τ₀)_+ ≤ φ₀ + φ(0)(1 − ρ²)²/(6ρ³).

*Proof.* Let c = √(1−ρ²), s = Σ a_j ε_j and f(x) = E_{g₀}(ρg₀ + x − τ₀)_+ =
ρψ((τ₀ − x)/ρ). f is C^∞ with f⁗(x) = ρ^{−3} φ″((τ₀ − x)/ρ), and
|φ″(u)| = |u² − 1| φ(u) ≤ φ(0) for all u: for |u| ≤ 1, (1 − u²)e^{−u²/2} ≤ 1; for
|u| > 1, (u² − 1)e^{−u²/2} is maximal at u² = 3 with value 2e^{−3/2} ≤ 1 (C8). So
sup|f⁗| ≤ φ(0)/ρ³. Now E(Z − τ₀)_+ = E f(cs). Replace ε_1, …, ε_r by independent
standard normals g_1, …, g_r one at a time. At step j, with R the (independent) sum of
the other terms, Taylor's formula to third order with fourth-order remainder gives
f(R + bξ) = Σ_{q=0}^{3} f^{(q)}(R) b^q ξ^q/q! + E_4, |E_4| ≤ sup|f⁗| b⁴ξ⁴/24, with
b = c a_j. Since ξ ∈ {ε_j, g_j} has Eξ = Eξ³ = 0 and Eξ² = 1 in both cases, the
expectations differ by at most sup|f⁗| b⁴ (Eε⁴ + Eg⁴)/24 = sup|f⁗| b⁴/6. Summing,
|E f(cs) − E f(c g′)| ≤ φ(0)c⁴ Σ_j a_j⁴/(6ρ³) ≤ φ(0)c⁴/(6ρ³), with g′ standard normal
(Σ a_j⁴ ≤ (Σ a_j²)² = 1). Finally E f(cg′) = E(ρg₀ + cg′ − τ₀)_+ = ψ(τ₀) because
ρg₀ + cg′ is standard normal. Add δ₀τ₀ and use Lemma 4.1. ∎

**Proposition 4.4 (Case A bound).** If α′ ≥ 3/100, then for every r, a, ρ there is
t ∈ {t_H, τ₀} ⊂ [0, t*] with Γ(ρ, a, t; α′) ≤ √α′ + 2φ₀.

*Proof.* Two certified facts (§8): (C3) 2D ≤ (1 − ρ₀)√(3/100), and (C4)
2b₀ ≤ √(3/100) where b₀ = φ(0)(1 − ρ₀)(1 + ρ₀)²/(6ρ₀³). Both right sides increase in
α′, so both hold with √α′ in place of √(3/100).

*ρ ≤ ρ₀.* With t = t_H, Lemma 4.2 gives Γ ≤ ρ√α′ + 2φ_H = ρ√α′ + 2φ₀ + 2D ≤
ρ√α′ + 2φ₀ + (1 − ρ₀)√α′ ≤ √α′ + 2φ₀.

*ρ ≥ ρ₀.* With t = τ₀, Lemma 4.3 gives Γ ≤ ρ√α′ + 2φ₀ + (1 − ρ)·2q(ρ)φ(0)/6, where
q(ρ) = (1 − ρ)(1 + ρ)²/ρ³ (we used (1−ρ²)² = (1−ρ)²(1+ρ)²). On (0, 1),
(ln q)′ = −1/(1−ρ) + 2/(1+ρ) − 3/ρ < 2 − 3 < 0, so q(ρ) ≤ q(ρ₀) and
2q(ρ)φ(0)/6 ≤ 2b₀ ≤ √α′. Hence Γ ≤ ρ√α′ + 2φ₀ + (1 − ρ)√α′ = √α′ + 2φ₀. ∎

**Corollary 4.5 (Case B bound).** For every α′ > 0 and all r, a, ρ:
Γ(ρ, a, t_H; α′) ≤ √α′ + 2φ_H. (Lemma 4.2 and ρ ≤ 1.)

## 5. The upper bound on E F

Let N ≥ 2 and M ≥ 1, and let X be an M × N matrix of iid uniform signs.

**Theorem 5.1.** Let b ∈ ℝ satisfy

  (H_b) for every r ∈ {1, …, M}, a ∈ S^{r−1} and ρ ∈ [0, 1] there is t ∈ [0, t*] with
  Γ(ρ, a, t; M/N) ≤ b.

Then for all ε ∈ (0, 1), η ∈ (0, 1] and β > 0, with r = ⌈εM⌉, J = binom(M, r),
𝓛 = (1 + 2/η)^r (π/(2η) + 2),

  E F ≤ N b + Err,  Err = Σ_{i=1}^{7} E_i,

  E₁ = 16√N + 4√(N ln(4J))                         (Step 2)
  E₂ = η√N ((21/5)√(M+N) + 2)                       (Step 3, net)
  E₃ = β^{−1} (N ln 2 + M ln(1 + 2/η))              (Step 3, soft maximum)
  E₄ = β² N M r^{−3/2}                              (Step 3, Lindeberg)
  E₅ = η√N (√M + √N + 2((21/5)√(r+N) + 2))          (Step 5, discretization)
  E₆ = 18√N + 8√(N ln(5𝓛))                          (Step 5, union bound)
  E₇ = 2t*                                          (Step 6)

*Proof.* Steps 1–6.

**Step 1 (dissection; deterministic).** For T ⊂ [M] let
B_T = {y ∈ ℝ^M : ‖y‖ ≤ 1, |y_μ| ≤ r^{−1/2} for μ ∉ T} and
F_T(A) = max_{σ∈Σ_k} sup_{y∈B_T} yᵀAσ for A ∈ ℝ^{M×N}. Then
max_{σ∈Σ_k} ‖Aσ‖ = max_{|T|=r} F_T(A). Indeed "≥" holds because B_T lies in the unit
ball. For "≤", take σ with u = Aσ ≠ 0 (if Aσ = 0 for all σ there is nothing to prove)
and let T index r largest |u_μ| (1 ≤ r ≤ M). If μ ∉ T, the r + 1 numbers u_ν², ν ∈ T ∪ {μ},
are each ≥ u_μ² and sum to at most ‖u‖², so |u_μ| ≤ ‖u‖/√(r+1); hence y = u/‖u‖ ∈ B_T
and yᵀAσ = ‖u‖.

**Step 2 (concentration of each piece and the union over T).** F_T is a supremum of
the linear functions A ↦ ⟨yσᵀ, A⟩_F with ‖yσᵀ‖_F = ‖y‖√N ≤ √N, so it is convex and
√N-Lipschitz in the Frobenius norm. Permuting rows preserves the law of X and maps F_T
to F_{πT}, so E F_T(X) = E F_{T₀}(X) for every |T| = r, T₀ := {1, …, r}. Let m_T be a
median of F_T(X). By T1, W_T := F_T(X) − m_T satisfies P(W_T ≥ s) ≤ 4e^{−s²/(16N)}, and
by Lemma 3.1, m_T ≤ E F_{T₀}(X) + 8√(πN). Lemma 3.2 with C = 4, v = 16N gives

  E F = E max_T F_T ≤ max_T m_T + E max_T (W_T)_+
      ≤ E F_{T₀}(X) + 8√(πN) + 4√(N ln(4J)) + 2√N/√(ln 4) ≤ E F_{T₀}(X) + E₁

(C8: 8√π + 2/√(ln 4) ≤ 16).

**Step 3 (Lindeberg on the capped rows only).** Let H be the matrix with the rows in
T₀ equal to those of X and the rows outside T₀ replaced by independent standard
Gaussian rows. Let Y ⊂ B_{T₀} be an η-net of B_{T₀} from Lemma 3.3, |Y| ≤ (1 + 2/η)^M,
K = Σ_k × Y, and f_β(A) = β^{−1} ln Σ_{(σ,y)∈K} exp(β yᵀAσ). For every A:

- (net) max_K yᵀAσ ≤ F_{T₀}(A) ≤ max_K yᵀAσ + η‖A‖_op √N, since for y ∈ B_{T₀} and
  y′ ∈ Y with ‖y − y′‖ ≤ η, (y − y′)ᵀAσ ≤ η ‖A‖_op ‖σ‖;
- (soft maximum) max_K ≤ f_β(A) ≤ max_K + β^{−1} ln|K|, and
  ln|K| ≤ ln binom(N, k) + M ln(1 + 2/η) ≤ N ln 2 + M ln(1 + 2/η);
- (third derivative) for an entry e = (μ, i) with μ ∉ T₀, writing ⟨·⟩ for the Gibbs
  average on K with weights exp(β yᵀAσ) and Z_e = y_μσ_i:
  ∂_e f_β = ⟨Z_e⟩, ∂_e² f_β = β(⟨Z_e²⟩ − ⟨Z_e⟩²), ∂_e³ f_β = β² ⟨(Z_e − ⟨Z_e⟩)³⟩.
  Since y ∈ B_{T₀}, |Z_e| ≤ r^{−1/2}, so |Z_e − ⟨Z_e⟩| ≤ 2r^{−1/2} and
  |∂_e³ f_β| ≤ 2r^{−1/2} β² ⟨(Z_e − ⟨Z_e⟩)²⟩ ≤ 2β² r^{−3/2}, for every value of all other
  entries. Entries in rows of T₀ are never differentiated.

Swap the N(M − r) entries outside T₀ from signs to Gaussians one at a time. At each
swap, the other entries form an independent family; expanding f_β in the swapped entry
x to second order with third-order remainder, |R(x)| ≤ 2β²r^{−3/2}|x|³/6, and using
that a sign and a standard Gaussian have the same first two moments, each swap changes
E f_β by at most (β² r^{−3/2}/3)(E|ε|³ + E|g|³) = β² r^{−3/2}(1 + 2√(2/π))/3 ≤ β² r^{−3/2}
(C8). (All terms are integrable: |f_β(A)| ≤ ‖A‖_F√N + ln|K|/β.) Hence
|E f_β(X) − E f_β(H)| ≤ β² N M r^{−3/2} = E₄. Chaining,

  E F_{T₀}(X) ≤ E max_K(X) + η√N E‖X‖_op ≤ E f_β(X) + E₂ ≤ E f_β(H) + E₂ + E₄
            ≤ E max_K(H) + E₂ + E₃ + E₄ ≤ E F_{T₀}(H) + E₂ + E₃ + E₄,

using Lemma 3.4 for E‖X‖_op ≤ (21/5)√(M+N) + 2.

**Step 4 (Sudakov–Fernique conditionally on the uncapped rows).** Write T = T₀,
y = (y_T, y′) with y′ ∈ ℝ^{M−r}, G for the Gaussian rows of H, and fix the value x of
X_T (the r × N sign rows); G is independent of X_T. On the index set Σ_k × B_T put

  𝒳(σ, y) = m(σ, y) + y′ᵀGσ,   𝒴(σ, y) = m(σ, y) + √N g·y′ + ‖y′‖ h·σ,

with common mean m(σ, y) = y_Tᵀ x σ, and g ∈ ℝ^{M−r}, h ∈ ℝ^N independent standard
Gaussian vectors. Then F_T(H) = sup 𝒳. With ρ = ‖y′‖, ρ̃ = ‖ỹ′‖:

  E(𝒳 − 𝒳̃)² = ‖y′σᵀ − ỹ′σ̃ᵀ‖_F² = N(ρ² + ρ̃²) − 2(y′·ỹ′)(σ·σ̃),
  E(𝒴 − 𝒴̃)² = N‖y′ − ỹ′‖² + ‖ρσ − ρ̃σ̃‖² = 2N(ρ² + ρ̃²) − 2N y′·ỹ′ − 2ρρ̃ σ·σ̃,

and the difference is

  E(𝒴 − 𝒴̃)² − E(𝒳 − 𝒳̃)² = N(ρ − ρ̃)² + 2(ρρ̃ − y′·ỹ′)(N − σ·σ̃) ≥ 0,

because y′·ỹ′ ≤ ρρ̃ and σ·σ̃ ≤ N. (Check: N(ρ−ρ̃)² + 2Nρρ̃ − 2Ny′·ỹ′ − 2ρρ̃σ·σ̃ +
2(y′·ỹ′)σ·σ̃ = N(ρ²+ρ̃²) − 2Ny′·ỹ′ − 2ρρ̃σ·σ̃ + 2(y′·ỹ′)σ·σ̃, which is the difference of the
two displayed lines.) Both processes are linear in y, hence continuous on the compact
B_T, and integrable. Corollary 3.6 (with S₀ = Σ_k, U = B_T) gives
E[F_T(H) | X_T = x] ≤ E sup 𝒴. To bound sup 𝒴, drop the caps (B_T lies in the unit
ball) and maximize over y with ‖y′‖ = ρ fixed: √N g·y′ ≤ ρ√N‖g‖ and
y_Tᵀxσ ≤ √(1−ρ²)‖xσ‖. With ρ = cos θ,

  sup 𝒴 ≤ max_{σ∈Σ_k} max_{θ∈[0,π/2]} [ cos θ (√N‖g‖ + h·σ) + sin θ ‖xσ‖ ].

Averaging over X_T,

  E F_T(H) ≤ E max_{σ∈Σ_k} max_{θ∈[0,π/2]} [ cos θ (√N‖g‖ + h·σ) + sin θ ‖X_Tσ‖ ].    (5.1)

**Step 5 (nets on the r uncapped rows and on the angle).** Let 𝒜 be an η-net of
S^{r−1} (Lemma 3.3, |𝒜| ≤ (1+2/η)^r) and Θ = {jη : 0 ≤ j ≤ ⌊π/(2η)⌋} ∪ {π/2}, so
|Θ| ≤ π/(2η) + 2 and |𝒜||Θ| ≤ 𝓛. For v = X_Tσ and a ∈ 𝒜 within η of v/‖v‖,
aᵀv ≥ ‖v‖ − η‖v‖, so ‖X_Tσ‖ ≤ max_{a∈𝒜} aᵀX_Tσ + η‖X_T‖_op √N. For θ ∈ [0, π/2] and
θ_j ∈ Θ with |θ − θ_j| ≤ η, |cos θ − cos θ_j| and |sin θ − sin θ_j| are ≤ η, and
|√N‖g‖ + h·σ| + |aᵀX_Tσ| ≤ √N(‖g‖ + ‖h‖ + ‖X_T‖_op). Therefore the right side of (5.1)
is at most E max_{(a,θ)∈𝒜×Θ} Φ_{a,θ} + E[η√N(‖g‖ + ‖h‖ + 2‖X_T‖_op)], where

  Φ_{a,θ} = cos θ √N ‖g‖ + max_{σ∈Σ_k} (cos θ h + sin θ X_Tᵀa)·σ.

By E‖g‖ ≤ √M, E‖h‖ ≤ √N and Lemma 3.4 for the r × N matrix X_T, the second expectation
is at most E₅.

*Concentration of each Φ_{a,θ}.* For fixed X_T, (g, h) ↦ Φ_{a,θ} is √(2N)-Lipschitz
(|ΔΦ| ≤ √N‖Δg‖ + √N‖Δh‖), so by T2, P(Φ − ψ_{a,θ}(X_T) ≥ s | X_T) ≤ e^{−s²/(4N)}, where
ψ_{a,θ}(X_T) := E[Φ_{a,θ} | X_T]. The function X_T ↦ Φ_{a,θ} is, for fixed (g, h), a
maximum of affine functions with ‖aσᵀ‖_F = √N, hence convex and √N-Lipschitz; so
ψ_{a,θ} is convex and √N-Lipschitz, and by T1 and (3.1),
P(ψ_{a,θ} − Eψ_{a,θ} ≥ s + 8√(πN)) ≤ 4e^{−s²/(16N)}. Since E Φ = E ψ,

  P(Φ_{a,θ} − EΦ_{a,θ} ≥ 2s + 8√(πN)) ≤ e^{−s²/(4N)} + 4e^{−s²/(16N)} ≤ 5e^{−s²/(16N)}.

Apply Lemma 3.2 (C = 5, v = 16N) to W_{a,θ} = (Φ_{a,θ} − EΦ_{a,θ} − 8√(πN))/2:

  E max_{a,θ} Φ_{a,θ} ≤ max_{a,θ} EΦ_{a,θ} + 8√(πN) + 8√(N ln(5𝓛)) + 4√N/√(ln 5)
                     ≤ max_{a,θ} EΦ_{a,θ} + E₆   (C8: 8√π + 4/√(ln 5) ≤ 18).

**Step 6 (a fixed direction is an iid problem).** Fix a ∈ S^{r−1}, θ, and put ρ = cos θ.
The coordinates w_i = Σ_{μ∈T} a_μ X_μi of w = X_Tᵀa are iid (the columns of X_T are iid)
with the law of Σ_μ a_μ ε_μ, and independent of h, so Z_i := cos θ h_i + sin θ w_i are
iid copies of Z_{ρ,a}, which is symmetric. For σ ∈ Σ_k with minus set S and any real t,
Z·σ = Σ_i Z_i + 2Σ_{i∈S}(−Z_i) and Σ_{i∈S} c_i ≤ kt + Σ_{i=1}^N (c_i − t)_+ for any reals c_i;
so max_{σ∈Σ_k} Z·σ ≤ Σ_i Z_i + 2kt + 2Σ_i (−Z_i − t)_+. Taking expectations, using
symmetry, k ≤ δ₀N + 1, and restricting to t ∈ [0, t*]:

  E Φ_{a,θ} ≤ ρ√N E‖g‖ + 2N(δ₀t + E(Z_{ρ,a} − t)_+) + 2t ≤ N Γ(ρ, a, t; M/N) + 2t*,

since ρ√N E‖g‖ ≤ ρ√N √(M − r) ≤ ρ√(NM). Choosing t by (H_b), E Φ_{a,θ} ≤ Nb + E₇.

Chaining Steps 2, 3, (5.1), 5 and 6 gives E F ≤ Nb + Σ E_i. ∎

**Entropy versus concentration (why the parameters are arranged this way).** Three
finite families appear. (i) The r-subsets T: ln J ≤ r ln(eM/r) ≤ (εM + 1) ln(e/ε)
is linear in N, but it is paid at the concentration exponent s²/(16N), so the cost
is s ≈ 4√(N ln J) ≈ 4N√(αε ln(e/ε)), which is o(N) as ε → 0. (ii) The net Y of
B_{T₀} has ln|Y| ≈ M ln(1 + 2/η); paying this at the concentration exponent would cost
≈ 4N√(α ln(1+2/η)), which does **not** vanish. So Y is never union-bounded: it only
enters the soft maximum, at cost ln|K|/β, and the Lindeberg error β²NMr^{−3/2} is
O(√N) for fixed (ε, β). (iii) The net 𝒜 × Θ has ln 𝓛 ≈ εM ln(1 + 2/η) and is
union-bounded at cost ≈ 8N√(αε ln(1+2/η)), which is o(N) as ε → 0 for fixed η.

## 6. The limit

**Corollary 6.1.** Let M = M(N) with M/N → α ∈ (0, ∞), and suppose b_N satisfies (H_{b_N})
for each N with b_N → b_∞. Then for every γ > 0 there is N_γ with
E F ≤ N(b_∞ + γ) for N ≥ N_γ.

*Proof.* For fixed (ε, η, β), as N → ∞: ln(4J)/N ≤ (ln 4 + (εM + 1) ln(e/ε))/N →
αε ln(e/ε); ln(5𝓛)/N → αε ln(1 + 2/η); r ≥ εM, so E₄/N ≤ β² ε^{−3/2} M^{−1/2} → 0;
E₇/N → 0. Hence

  limsup_N (E F/N − b_N) ≤ R(ε, η, β) := 4√(αε ln(e/ε)) + 8√(αε ln(1+2/η))
      + η(√α + 1 + (63/5)√(1+α)) + β^{−1}(ln 2 + α ln(1 + 2/η)).

(The η term collects E₂/N → (21/5)η√(1+α) and E₅/N ≤ η(√(M/N) + 1 + (42/5)√((r+N)/N) +
4/√N) → at most η(√α + 1 + (42/5)√(1+α)).) The left side does not depend on
(ε, η, β), and inf R = 0: given γ, choose η with η(√α + 1 + (63/5)√(1+α)) < γ/3, then,
for that η, ε so small that the two ε-terms are < γ/3, and β so large that the β-term is
< γ/3. So limsup (E F/N − b_N) ≤ 0, and b_N → b_∞ gives the claim. ∎

The order of limits is therefore: N → ∞ first (this is what kills the Lindeberg term
β²ε^{−3/2}M^{−1/2}, and the √N terms); then β → ∞ and ε → 0 at fixed η (in either
order); then η → 0.

**Corollary 6.2 (hypothesis (U) of PROOF.md, at τ₀).** If M/N → α > 3/100 then
E F ≤ √(NM) + 2Nφ₀ + o(N). If M/N → α > 0, then E F ≤ √(NM) + 2Nφ_H + o(N).

*Proof.* Case A: for N large, α_N := M/N ≥ 3/100, so by Proposition 4.4, (H_b) holds
with b_N = √α_N + 2φ₀ → √α + 2φ₀. Case B: Corollary 4.5 gives (H_b) with
b_N = √α_N + 2φ_H for every N. Apply Corollary 6.1. ∎

## 7. Proof of Theorem 2′

Fix ν; X = X^{(ν)} has M = P − 1 rows and M/N → α.

**The bound b.** In Case A (α ∈ (3/100, 0.12948)) let b_∞ = √α + 2φ₀; in Case B
(α ∈ (0, 3/100]) let b_∞ = √α + 2φ_H. In both cases

  s_α := α + 4δ₀(1 − δ₀) − b_∞² > 0.

Indeed, for p > 0 with δ₀(1−δ₀) > p², (√α + 2p)² < α + 4δ₀(1−δ₀) ⟺ 4p√α + 4p² < 4δ₀(1−δ₀)
⟺ √α < (δ₀(1−δ₀) − p²)/p. With p = φ₀ the right side is √A(τ₀), and A(τ₀) > 0.12948
(C5). With p = φ_H = δ₀√(2L₀) the right side is (1 − δ₀ − 2δ₀L₀)/√(2L₀) = √A_H, and
A_H > 3/100 (C6); the numerator is certified positive.

**Upper tail of F.** F = max_{σ∈Σ_k} ‖Xσ‖ is convex and √N-Lipschitz in X
(|F(X) − F(X′)| ≤ max_σ ‖(X − X′)σ‖ ≤ ‖X − X′‖_F √N). By (3.1),
P(F ≥ E F + 8√(πN) + s) ≤ 4e^{−s²/(16N)}. Fix γ > 0 (chosen below). By Corollary 6.2,
E F ≤ N(b_∞ + γ) for N ≥ N_γ; for N also so large that 8√(πN) ≤ γN, taking s = γN,

  P(F ≥ N(b_∞ + 3γ)) ≤ 4e^{−γ²N/16}.

**Lower tail of ‖X**1**‖².** ‖X**1**‖² = Σ_μ S_μ², where S_μ = Σ_i X_μi are independent,
E S_μ² = N and E S_μ⁴ = 3N² − 2N ≤ 3N². For ε′ ∈ (0, 1) and λ > 0, by Markov's
inequality and e^{−x} ≤ 1 − x + x²/2 (x ≥ 0),

  P(‖X**1**‖² ≤ (1−ε′)MN) ≤ e^{λ(1−ε′)MN} (E e^{−λS_1²})^M,
  E e^{−λS_1²} ≤ 1 − λN + (3/2)λ²N² ≤ exp(−λN + (3/2)λ²N²).

With λN = ε′/3 the bound is exp(M(−ε′λN + (3/2)λ²N²)) = exp(−ε′²M/6).

**Combination.** For N ≥ 1/(1 − 2δ₀), k/N ∈ [δ₀, δ₀ + 1/N] ⊂ [δ₀, 1 − δ₀], and x ↦ x(1−x)
is at least δ₀(1−δ₀) there, so 4k(N−k) ≥ 4δ₀(1−δ₀)N². On the event
{F < N(b_∞ + 3γ)} ∩ {‖X**1**‖² > (1−ε′)MN}, (★) holds as soon as

  (b_∞ + 3γ)² ≤ (1 − ε′) M/N + 4δ₀(1 − δ₀).

Choose γ > 0 and ε′ ∈ (0, 1) with (b_∞ + 3γ)² − b_∞² ≤ s_α/3 and ε′(α + 1) ≤ s_α/3. For N
so large that |M/N − α| ≤ s_α/6, the right side is at least α + 4δ₀(1−δ₀) − ε′α −
(1−ε′)s_α/6 ≥ b_∞² + s_α − s_α/3 − s_α/6 > (b_∞ + 3γ)². So for N ≥ N₀ (N₀ depending on
α, the sequence M(N), γ and ε′, but not on ν),

  P((★) fails for ν) ≤ 4e^{−γ²N/16} + e^{−ε′²M/6} ≤ 5e^{−cN},
  c = min(γ²/16, ε′²α/12)   (using M ≥ αN/2 for large N).

**Union over patterns.** The law of X^{(ν)} is the same for every ν, so the same N₀ and c
serve all ν, and P(some ν fails) ≤ 5P e^{−cN}. On the complement, §2 gives statements 1
and 2 for every pattern simultaneously. Case A and Case B cover (0, 0.12948). ∎

**Remarks.**

- Only the constants of §4 depend on the patterns being signs (Hoeffding for the r
  uncapped rows, and the fourth moment Eε⁴ = 1 ≤ 3). Talagrand's inequality uses that
  the entries are bounded; the Lindeberg steps use the first two (Step 3) or three
  (Lemma 4.3) moments.
- The same proof with the Hoeffding route alone at δ_H = Q(2.2217) gives the barrier for
  α < A_H(δ_H) = 0.0879656… (C7). This is the best one can do without Lemma 4.3, but it
  is not needed: at the radius δ₀ the Hoeffding route already reaches 0.0801 > 0.03.
- In model R the target pattern is one of the P stored patterns, so Theorem 2′ is a
  statement about every stored pattern, not one target with independent interference.

## 8. Certified numbers

`check_pm1.py` (same folder) uses `fractions.Fraction` only in the certified path:
π by Machin's formula, e^{−x} and the normal integral by Taylor series with explicit
tail bounds, ln by the atanh series with an explicit tail bound, square roots by
integer-square-root brackets. Every inequality below is decided on rational interval
endpoints. A python-flint (arb, 256 bits) cross-check of the main quantities runs
alongside and must agree; it is not part of the decision. Output (decimals rounded
outward):

| check | claim | certified |
| :-- | :-- | :-- |
| C1 | 0 < δ₀ < 1/2; L₀ > 1/2; τ₀² ≤ 2L₀ | δ₀ ∈ [0.033932663104, 0.033932663105], L₀ ∈ [3.383377215072, 3.383377215073] |
| C2 | D > 0 | D ∈ [0.012940278901, 0.012940278902] |
| C3 | (a): (2D/(1−ρ₀))² ≤ 3/100 | threshold ≤ 0.029769034319 |
| C4 | (b): (2b₀)² ≤ 3/100 | b₀ ≤ 0.055582330742, threshold ≤ 0.012357581963 |
| C5 | A(τ₀) > 0.12948 | A(τ₀) ∈ [0.129489982621, 0.129489982622] |
| C6 | A_H(δ₀) > 3/100 | A_H ∈ [0.080151204746, 0.080151204747] |
| C7 | A_H(Q(22217/10000)) > 0.0879 (informational) | [0.087965649827, 0.087965649828] |
| C8 | 2e^{−3/2} ≤ 1; 2√(2 ln 9) ≤ 21/5; (1 + 2√(2/π))/3 ≤ 1; 8√π + 2/√(ln 4) ≤ 16; 8√π + 4/√(ln 5) ≤ 18; 1 − 2δ₀ ≥ 0.932 | all pass |
| C9 | consistency: barrier inequality at α = 0.12948 (Case A) and α = 0.03 (Case B) | margins ≥ 4.1·10⁻⁶ and ≥ 0.0388 |

(C3 and C4 are checked at α = 3/100 and also at 0.12948; by monotonicity they hold on
[3/100, ∞).) Floating-point sanity checks (not certified, not used): for a single sign
(r = 1) and α ∈ {0.005, 0.03, 0.1294}, a grid over ρ of ρ√α + 2 inf_t(δ₀t + E(Z−t)_+)
never exceeds √α + 2φ₀ (maximum at ρ = 1, excess ≤ 6·10⁻¹⁷) and stays below √α + 2φ_H by
0.0259; this agrees with Proposition 4.4 and Corollary 4.5.

## 9. Changes from the draft (`pm1-plan.md`)

I found no mathematical error in the draft's Steps 1–6 or in its Proposition 4. The
changes are a gap closed, statements narrowed to what is used, and steps that the draft
cited or sketched written out.

1. **Small α (the gap).** The draft covered only α ≥ 0.03 and, below that, pointed to an
   uncertified grid search over τ or to classical bounds. Following the caller's fix, a
   Hoeffding-only bound covers small α. Change from that fix: I apply it at the **same**
   radius δ₀ = Q(τ₀) instead of at τ_H ≈ 2.2217. A_H(δ₀) = 0.08015 > 0.03 (C6), so one
   radius ⌈δ₀N⌉ and one overlap bound 0.932 serve every α ∈ (0, 0.12948). The τ_H version
   (0.087966, C7) is certified but not needed. In the Hoeffding route τ plays no role;
   it is parametrized by δ alone, A_H(δ) = (1 − δ − 2δ ln(1/δ))²/(2 ln(1/δ)).
2. **Theorem 3 restated without V and Ψ_δ.** The draft bounded E F by V(α, δ), a sup of an
   inf over t, and asserted without proof that the minimizing t is at most √(2 ln(1/δ)).
   (That claim is true, because P(Z > √(2 ln(1/δ))) ≤ δ by Hoeffding, but it is not
   needed.) Theorem 5.1 instead takes the hypothesis (H_b): one t ∈ [0, t*] per (ρ, a).
   Proposition 4.4 supplies t ∈ {t_H, τ₀} explicitly, and Lemma 4.1 replaces the
   minimization at ρ = 1 by a direct evaluation.
3. **Step 2.** A median per T (equal laws make them interchangeable, but this avoids the
   point); mean-versus-median proved (Lemma 3.1); the union-bound integration proved
   (Lemma 3.2). Constant 16√N + 4√(N ln 4J) (draft: 17√N + 4√(N ln 4J)).
4. **Step 3.** The draft's E‖X‖_op ≤ √(π/2)(√N + √M) (contraction plus Gordon) is correct
   but cited; replaced by the elementary Lemma 3.4, (21/5)√(M+N) + 2. Third-moment constant
   explicit: (1 + 2√(2/π))/3 ≤ 1. The draft's parenthetical appeal to Han–Shen (2023,
   Thm 2.5) is not used.
5. **Step 4.** Sudakov–Fernique with a common non-zero mean is proved in full by
   interpolation (Theorem 3.5) and extended to the compact index set (Corollary 3.6). The
   increment identity was rechecked term by term.
6. **Step 5.** The draft's "Gaussian concentration and Talagrand" is made explicit as a
   two-stage bound (Gaussian in (g, h) given X_T; Talagrand for the conditional mean, which
   is convex and √N-Lipschitz), with constant 5e^{−s²/(16N)}; the discretization error is
   written out (E₅).
7. **Limit order.** Made precise: N → ∞ innermost, then β → ∞ and ε → 0 at fixed η, then
   η → 0 (Corollary 6.1). Added the entropy-versus-exponent accounting at the end of §5,
   including why the net Y must go through the soft maximum and never through a union
   bound.
8. **From expectation to every pattern.** The draft's "Consequence" relied on the sketch in
   PROOF.md §6. It is written out here: concentration of F with explicit exponent, the
   Chernoff lower tail of ‖X**1**‖² with the moment computation, the choice of γ and ε′
   from the slack s_α, and the union over all P patterns with an N₀ independent of ν.
9. **Overstated remark removed.** The draft said the proof "covers any symmetric entries
   with variance 1 and variance proxy 1." As written it uses bounded entries (Talagrand, in
   Step 2, Step 5 and the upper tail of F) and Eε⁴ ≤ 3 (Lemma 4.3), so the remark needs at
   least a replacement for T1; it is not claimed here.
10. **Out of scope.** The draft's §6 (reference-state extension to 0.1347) and its
    Proposition 6 were not reviewed and are not part of this proof.

## 10. What would most likely hide an error

For a second reader, in order: (1) the increment identity of Step 4 and the use of
Corollary 3.6 conditionally on X_T; (2) Step 3's third-derivative bound, which uses the
cap only for swapped rows; (3) the claim in Step 6 that X_Tᵀa has iid coordinates (it
does: the columns of X_T are iid); (4) the sign bookkeeping in Lemma 4.3 and
sup|φ″| = φ(0); (5) that the cited forms of T1 and T2 are the ones in the books.

## 11. References

- M. Talagrand (1988), An isoperimetric theorem on the cube and the Khintchine–Kahane
  inequalities, Proc. AMS 104, 905–909. M. Ledoux (2001), The Concentration of Measure
  Phenomenon, AMS, Corollary 4.10.
- S. Boucheron, G. Lugosi, P. Massart (2013), Concentration Inequalities, OUP,
  Theorem 5.6.
- X. Fernique (1975), LNM 480. R. A. Vitale (2000), Proc. AMS 128, 3043–3046.
  S. Chatterjee (2005), An error bound in the Sudakov–Fernique inequality,
  arXiv:math/0510424.
- S. Oymak, J. A. Tropp (2018), Universality laws for randomized dimension reduction,
  Inf. Inference 7 (arXiv:1511.09433), §9: the dissection device of Step 1.
- R. Vershynin (2018), High-Dimensional Probability, CUP, §4.2 (nets), §4.4 (operator
  norm by nets).
- J. Feng, M. Shcherbina, B. Tirozzi (2001), On the critical capacity of the Hopfield
  model, Comm. Math. Phys. 216 (the 0.113 bound).

Theorem and section numbers were taken from memory and were not re-checked against the
books in this session.
