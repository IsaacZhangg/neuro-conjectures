# Reference-state refinement of Theorem 1 (model G)

This file proves Theorem 3, which raises the certified load for model G
(one ±1 target, independent N(0,1) interference patterns) from
0.12948 (Theorem 1 of [PROOF.md](PROOF.md)) to
**0.134716**. The only change of method is the reference point of the energy
barrier: instead of comparing the Hamming sphere Σ_k with H(ξ), it is compared
with H(σ₁), where σ₁ is ξ with a small, explicitly chosen set of spins
flipped. Nothing here is claimed to be new mathematics.

Status: every step of Theorem 3 is written out below. Three ingredients are
used without re-proof: Lemmas 1, 3 and 4 of PROOF.md (proved there), and the
standard tail inequalities listed in Section 5 with citations. Theorem and
page numbers of the cited books were taken from memory and not re-checked in
this session (the same caveat as PROOF.md, Section 7). Section 8 lists what
is not proved.

## 1. Setting

Notation as in PROOF.md, Section 1. N ≥ 2 spins, target ξ ∈ {±1}^N, M ≥ 0
interference patterns forming the M × N matrix G with independent N(0, 1)
entries, J_ij = (1/N)(ξ_iξ_j + Σ_μ G_μiG_μj) for i ≠ j, J_ii = 0,
h_i(σ) = Σ_{j≠i} J_ij σ_j, H(σ) = −½ Σ_{i≠j} J_ij σ_iσ_j. Stable fixed point:
σ_i h_i(σ) > 0 for all i. φ, Φ are the standard normal density and
distribution function, Q = 1 − Φ. **1** is the all-ones vector, e_T ∈ {0,1}^N
the indicator vector of T ⊂ [N], and P = I − **1** **1**ᵀ/N. Σ_k is the set of
σ with exactly k coordinates equal to −1 (after the gauge ξ = **1**).

Used from PROOF.md:

- **Lemma 1.** H(σ) = −(1/2N)[(ξ·σ)² + ‖Gσ‖²] + (1/2N)[N + Σ_{μ,i} G_μi²].
  Hence H(σ) > H(ρ) ⟺ (ξ·σ)² + ‖Gσ‖² < (ξ·ρ)² + ‖Gρ‖².
- **Lemma 3.** In model G (N ≥ 2), almost surely no two configurations at
  Hamming distance 1 have equal energy.
- **Lemma 4.** For τ > 0, δ = Q(τ), k = ⌈δN⌉ and F := max_{σ∈Σ_k} ‖Gσ‖,
  E F ≤ √(NM) + 2Nφ(τ) + 2τ.
- **Gauge.** σ_i ↦ ξ_iσ_i, G_μi ↦ G_μiξ_i preserves the law of G, energies and
  Hamming distances, and maps ξ to **1**. From Section 3 on, ξ = **1**.

## 2. The constant

For 0 < τ < θ put δ = Q(τ) and

  D(τ, θ) = φ(τ) − φ(θ)(1 − 2Q(θ)),
  Num(τ, θ) = δ(1 − δ) − φ(τ)² + φ(θ)² − Q(θ)(1 − Q(θ)),
  A₂(τ, θ) = (Num/D)²  when Num > 0 (and 0 otherwise).

Since θ > τ > 0 we have 0 < φ(θ) < φ(τ) and 0 < 1 − 2Q(θ) < 1, so D > 0
automatically, and Q(θ) < Q(τ) = δ. Also, for m ≥ 0,

  (2.1)  Num − √m D > 0  ⟺  √m < Num/D  ⟺  m < A₂  (when Num > 0).

As θ → ∞, Num → δ(1−δ) − φ(τ)² and D → φ(τ), so A₂(τ, θ) → A(τ), the
constant of Theorem 1; Theorem 3 below therefore contains Theorem 1 in the
limit, and sup A₂ ≥ sup A.

## 3. Generalised barrier lemma

**Lemma 2\* (barrier relative to an interior reference point).** Fix k ≥ 1
and a configuration ρ with d(ρ, ξ) ≤ k − 1. Suppose that H(σ) > H(ρ) for every
σ with d(σ, ξ) = k, and that H(σ) ≠ H(σ′) whenever σ, σ′ differ in exactly
one spin. Then there is a stable fixed point σ* with d(σ*, ξ) ≤ k − 1 and
H(σ*) ≤ H(ρ).

*Proof.* B = {σ : d(σ, ξ) ≤ k − 1} is finite and contains ρ. Let σ* minimise
H on B; then H(σ*) ≤ H(ρ). Fix i and let σ′ be σ* with spin i flipped, so
d(σ′, ξ) ≤ k. If d(σ′, ξ) ≤ k − 1, then σ′ ∈ B, H(σ′) ≥ H(σ*) by minimality
and H(σ′) ≠ H(σ*) by the no-ties assumption, so H(σ′) > H(σ*). If
d(σ′, ξ) = k, then H(σ′) > H(ρ) ≥ H(σ*) by hypothesis. In both cases
H(σ′) − H(σ*) > 0. The terms of H containing σ_i sum to −σ_i h_i(σ) (h_i does
not depend on σ_i because J_ii = 0), so flipping spin i changes H by exactly
2σ*_i h_i(σ*). Hence σ*_i h_i(σ*) > 0 for every i. ∎

With ρ = ξ this is Lemma 2 of PROOF.md. As there, dropping the no-ties
assumption gives σ*_i h_i(σ*) ≥ 0 for all i.

## 4. Gaussian decomposition and the exact energy identity

Assume ξ = **1** and M ≥ 1. Put s = G**1** ∈ ℝ^M and Y = GP, so
G = s**1**ᵀ/N + Y.

**Lemma 5 (decomposition).** Let ν ~ N(0, 1/N) be independent of G
(auxiliary randomisation; it changes no statement about G, see Step 0 of the
proof of Theorem 3).

(a) s ~ N(0, N I_M) and Y is independent of s.

(b) For every fixed unit vector u ∈ ℝ^M, the vector uᵀY ∈ ℝ^N (a row vector,
identified with its transpose) has law N(0, P), and it is independent of
(I − uuᵀ)Y. For every fixed e ∈ ℝ^N, (I − uuᵀ)Ye ~ N(0, ‖Pe‖²(I − uuᵀ)).

(c) On {s ≠ 0} (probability one) define ŝ = s/‖s‖, z = Yᵀŝ ∈ ℝ^N,
R = (I − ŝŝᵀ)Y and w = z + ν**1**. Then, conditionally on s: w ~ N(0, I_N),
R is independent of w, and for each fixed e, Re ~ N(0, ‖Pe‖²(I − ŝŝᵀ)).
Moreover Y = ŝzᵀ + R, Pw = z, and (1/N)**1**ᵀw = ν.

*Proof.* (a) The rows g_μ of G are independent N(0, I_N). For each μ,
(g_μ·**1**, Pg_μ) is jointly Gaussian with Cov(Pg_μ, g_μ·**1**) = P**1** = 0, so
g_μ·**1** and Pg_μ are independent; across μ everything is independent. So
s = (g_μ·**1**)_μ is independent of Y (whose rows are Pg_μ), and
s_μ ~ N(0, ‖**1**‖²) = N(0, N).

(b) Let H_u be an orthogonal M × M matrix with first row uᵀ (e.g. a Householder
reflection, which depends measurably on u). The columns of G are independent
N(0, I_M), so the columns of H_uG are independent N(0, H_uH_uᵀ) = N(0, I_M):
H_uG again has independent N(0, 1) entries, and its rows r_1, …, r_M are
independent N(0, I_N). Now uᵀY = r_1P, and
(I − uuᵀ)Y = H_uᵀ diag(0, 1, …, 1) H_u G P is a function of r_2, …, r_M only.
Hence uᵀY ⫫ (I − uuᵀ)Y, and uᵀY = r_1P ~ N(0, PPᵀ) = N(0, P). For fixed e,
GPe has independent N(0, ‖Pe‖²) coordinates (g_μ·Pe), so
(I − uuᵀ)GPe ~ N(0, ‖Pe‖²(I − uuᵀ)(I − uuᵀ)ᵀ) = N(0, ‖Pe‖²(I − uuᵀ)).

(c) By (a), conditioning on s = s₀ ≠ 0 leaves the law of (Y, ν) unchanged and
fixes u = s₀/‖s₀‖, so (b) applies with that u (formally: by independence the
regular conditional law of (Y, ν) given s is the unconditional one, and
(z, R) are fixed measurable functions of (s, Y)). Given s, ν is independent of
(z, R), and (z, ν) is independent of R, so w = z + ν**1** is independent of R.
Its law is N(0, P + **1****1**ᵀ/N) = N(0, I_N). Finally
ŝŝᵀY = ŝ(Yᵀŝ)ᵀ = ŝzᵀ, so Y = ŝzᵀ + R; **1**ᵀz = ŝᵀY**1** = ŝᵀGP**1** = 0, so
Pz = z, Pw = z + νP**1** = z, and **1**ᵀw/N = **1**ᵀz/N + ν = ν. ∎

**Definition (reference state).** Fix θ > 0. Let
T = {i : w_i < −θ}, σ₁ = **1** − 2e_T, p = |T|/N, and
Nζ := −Σ_{i∈T} z_i.

T is a measurable function of w, hence (by Lemma 5(c)) independent of R given s.
Note d(σ₁, **1**) = |T| = pN and **1**·σ₁ = N(1 − 2p).

**Lemma 6 (exact energy identity).** On {s ≠ 0},

  Gσ₁ = ŝ (‖s‖(1 − 2p) + 2Nζ) − 2Re_T,
  ‖Gσ₁‖² = (‖s‖(1 − 2p) + 2Nζ)² + 4‖Re_T‖²,
  (**1**·σ₁)² + ‖Gσ₁‖² = N²(1 − 2p)² + (‖s‖(1 − 2p) + 2Nζ)² + 4‖Re_T‖².

*Proof.* G**1** = s. G e_T = s(**1**·e_T)/N + Ye_T = ps + ŝ(z·e_T) + Re_T
= ps − Nζ ŝ + Re_T. So Gσ₁ = s − 2ps + 2Nζŝ − 2Re_T, and s = ‖s‖ŝ gives the
first line. Since ŝᵀR = ŝᵀ(I − ŝŝᵀ)Y = 0, Re_T ⊥ ŝ and Pythagoras gives the
second line; the third adds (**1**·σ₁)² = N²(1 − 2p)². ∎

(The identity was also checked numerically, relative error ≤ 6·10⁻¹⁶, by
`mc_refstate.py`.)

Two further exact facts used below:

- ‖Pe_T‖² = |T| − |T|²/N = N p(1 − p).
- Nζ = −Σ_{i∈T}(w_i − ν) = Σ_{i=1}^N Y_i + |T|ν, where
  Y_i := (−w_i) 1{−w_i > θ} (using z = Pw = w − ν**1**).

**Lemma 7 (conditional law of ‖Re_T‖²).** Conditionally on (s, w), on the
event {0 < |T| < N}, X := ‖Re_T‖²/‖Pe_T‖² has the χ² law with M − 1 degrees
of freedom (X = 0 when M = 1).

*Proof.* Given s, R ⫫ w (Lemma 5(c)), so the conditional law of R given
(s, w) is its conditional law given s. Since e_T is a function of w, given
(s, w) the vector Re_T is N(0, ‖Pe_T‖²(I − ŝŝᵀ)) by Lemma 5(c) applied to the
fixed vector e = e_T. I − ŝŝᵀ is an orthogonal projection of rank M − 1, so
‖Re_T‖²/‖Pe_T‖² ~ χ²_{M−1} when ‖Pe_T‖ > 0, i.e. when 0 < |T| < N. ∎

## 5. Tail inequalities used

All are for fixed N, M; "χ²_D" denotes a chi-square variable with D ≥ 1
degrees of freedom.

- **(T1) Laurent–Massart** (Ann. Statist. 28 (2000), Lemma 1):
  P(χ²_D ≤ D − 2√(Dx)) ≤ e^{−x} for x > 0.
- **(T2) Gaussian Lipschitz concentration** (Tsirelson–Ibragimov–Sudakov;
  Boucheron–Lugosi–Massart, *Concentration Inequalities*, Thm 5.6): if f is
  L-Lipschitz on ℝ^n (Euclidean norm) and g ~ N(0, I_n), then
  P(f(g) ≥ E f(g) + r) ≤ e^{−r²/(2L²)}. Used only for F, which is
  √N-Lipschitz in G for the Frobenius norm (PROOF.md, proof of Theorem 1), so
  P(F ≥ EF + εN) ≤ e^{−ε²N/2}.
- **(T3) Hoeffding** (1963): for independent Bernoulli(q) variables,
  P(|(1/N)Σ B_i − q| > ε) ≤ 2e^{−2ε²N}.
- **(T4) Lower tail of a sum of nonnegative i.i.d. variables.** If Y_i ≥ 0
  are i.i.d. with E Y = μ and E Y² ≤ v, then for ε > 0,
  P(Σ_{i≤N} Y_i ≤ N(μ − ε)) ≤ exp(−Nε²/(2v)).
  *Proof (complete).* For x ≥ 0, e^{−x} ≤ 1 − x + x²/2 (the Taylor remainder
  −x³e^{−ξ}/6 is ≤ 0). So for λ ≥ 0, E e^{−λY} ≤ 1 − λμ + λ²v/2
  ≤ exp(−λμ + λ²v/2). By Markov's inequality applied to e^{−λΣY},
  P(ΣY ≤ N(μ − ε)) ≤ e^{λN(μ−ε)} (E e^{−λY})^N ≤ exp(N(−λε + λ²v/2)).
  Take λ = ε/v. ∎
  This is the tool for the sum Σ_{i∈T} z_i, which is *not* a Lipschitz
  function of z: by writing it as Σ_i Y_i + |T|ν with T defined through the
  i.i.d. vector w, it becomes an i.i.d. sum of nonnegative terms plus a
  scalar Gaussian term. Here Y_i = X_i 1{X_i > θ} with X_i = −w_i ~ N(0, 1)
  i.i.d., μ = ∫_θ^∞ xφ(x)dx = φ(θ), and E Y² = θφ(θ) + Q(θ) ≤ E X² = 1, so
  P(Σ Y_i ≤ N(φ(θ) − ε)) ≤ e^{−ε²N/2}.
- **(T5) Scalar Gaussian tail.** For ν ~ N(0, 1/N), P(|ν| > ε) ≤ 2e^{−ε²N/2}
  (Chernoff: P(ν > ε) ≤ e^{−λε} E e^{λν} = e^{−λε + λ²/(2N)}, λ = εN).

## 6. Theorem 3 and its proof

**Theorem 3.** In model G, fix 0 < τ < θ with Num(τ, θ) > 0, put δ = Q(τ) and
k = ⌈δN⌉, and let M = M(N) satisfy limsup M/N < A₂(τ, θ). There are c > 0 and
N₀ such that for every N ≥ N₀, with probability at least 1 − 8e^{−cN}:

1. there is a configuration ρ with d(ρ, ξ) ≤ k − 1 such that H(σ) > H(ρ) for
   every σ with d(σ, ξ) = k; and
2. there is a stable fixed point σ* with d(σ*, ξ) ≤ k − 1, hence
   ξ·σ*/N ≥ 1 − 2δ − 2/N.

**Corollary.** At τ₀ = 9711/5000 = 1.9422 and θ₀ = 1582/625 = 2.5312,
[check_refstate.py](check_refstate.py) certifies, in exact rational arithmetic,
Q(θ₀) ∈ [0.005683650576233, 0.005683650576234] < Q(τ₀) ∈
[0.026056446551134, 0.026056446551135], Num > 0, D > 0 and
A₂(τ₀, θ₀) ∈ [0.134716551489, 0.134716551490]. So for every load with
limsup M/N < 0.134716, with probability 1 − O(e^{−cN}) the target of model G
has a stable fixed point within distance ⌈0.02606N⌉ − 1 of it (overlap at
least 0.9478 − 2/N). Numerically (floating point, not certified, not needed)
sup_{τ<θ} A₂ ≈ 0.1347166 at τ ≈ 1.9422, θ ≈ 2.5312.

*Proof of Theorem 3.*

*Step 0 (reductions).* By the gauge, ξ = **1**. If M = 0, then
h_i(**1**) = (N − 1)/N > 0 for all i, so ξ is a stable fixed point and both
statements hold with ρ = σ* = ξ (k ≥ 1). Assume M ≥ 1 from now on. We work on
the product space of G and an independent ν ~ N(0, 1/N). The events in
statements 1 and 2 depend on G only; if an event E of the product space is
contained in (statement holds) × ℝ, then P_G(statement) ≥ P(E). By Lemma 3,
outside a null set of G there are no single-flip ties, and then statement 2
follows from statement 1 by Lemma 2\*. So it suffices to show
P(E) ≥ 1 − 8e^{−cN} for an event E on which statement 1 holds with ρ = σ₁.

Choose a with limsup M/N < a < A₂(τ, θ), and N₁ with M ≤ aN for N ≥ N₁; write
m = M/N ∈ (0, a]. By (2.1) the **slack**

  S := 4(Num − √a D) > 0,

and since D > 0, 4(Num − √m D) ≥ S for every m ∈ [0, a].

*Step 1 (reduction to one inequality).* For σ ∈ Σ_k, **1**·σ = N − 2k and
‖Gσ‖ ≤ F. By Lemma 1, H(σ) > H(σ₁) for all σ ∈ Σ_k as soon as

  (★★)  (N − 2k)² + F² < N²(1 − 2p)² + ‖Gσ₁‖²,

and σ₁ is admissible as ρ as soon as |T| ≤ k − 1.

*Step 2 (the good event).* Fix ε, η ∈ (0, 1) (chosen in Step 5) with

  (6.1)  ε < min(Q(θ), δ − Q(θ), φ(θ)/2, 1/2 − Q(θ)).

Let E be the intersection of

- E1 = {F ≤ EF + εN}, P(E1ᶜ) ≤ e^{−ε²N/2} by (T2);
- E2 = {‖s‖²/N ≥ M − 2√(MηN)}, P(E2ᶜ) ≤ e^{−ηN} by (T1) with x = ηN, since
  ‖s‖²/N ~ χ²_M (Lemma 5(a));
- E3 = {|p − Q(θ)| ≤ ε}, P(E3ᶜ) ≤ 2e^{−2ε²N} by (T3), since, given s, w is
  standard Gaussian (Lemma 5(c)), so the 1{w_i < −θ} are i.i.d. Bernoulli(Q(θ))
  (also unconditionally, as the conditional law does not depend on s);
- E4 = {Σ_i Y_i ≥ N(φ(θ) − ε)}, P(E4ᶜ) ≤ e^{−ε²N/2} by (T4);
- E5 = {|ν| ≤ ε}, P(E5ᶜ) ≤ 2e^{−ε²N/2} by (T5);
- E6 = {X ≥ (M − 1) − 2√((M − 1)ηN)} ∪ E3ᶜ (X as in Lemma 7). On E3, (6.1)
  gives 0 < |T| < N, so by Lemma 7 and (T1),
  P(E6ᶜ) = E[P(E3 ∩ {X < …} | s, w)] ≤ e^{−ηN} (for M = 1 the event holds
  surely: X = 0 and the bound is 0).

Hence P(Eᶜ) ≤ 2e^{−ηN} + 6e^{−min(ε²/2, η)N} ≤ 8e^{−cN} with
c = min(ε²/2, η).

*Step 3 (consequences on E, deterministic).*

(i) |T| ≤ N(Q(θ) + ε) < Nδ ≤ k, so |T| ≤ k − 1: σ₁ lies in the open ball.
Also p ≤ Q(θ) + ε < 1/2.

(ii) F/N ≤ √m + 2φ(τ) + ε₁ with ε₁ := ε + 2τ/N, by E1 and Lemma 4.

(iii) b := ‖s‖/N ≥ b₀ := √(max(0, m − 2√(mη))), by E2.

(iv) ζ ≥ c₀ := φ(θ) − 2ε > 0, because Nζ = Σ Y_i + |T|ν ≥ N(φ(θ) − ε) − N|ν|
by E4 and E5.

(v) X/N ≥ x₀ := m − 1/N − 2√(mη), by E6 (on E3) and (M − 1)/N ≤ m.

(vi) k/N ∈ [δ, δ + 1/N]. For N ≥ N₂ := ⌈1/(1 − 2δ)⌉ + 1 this lies in
[δ, 1 − δ], where x(1 − x) ≥ δ(1 − δ), so (N − 2k)²/N² = 1 − 4(k/N)(1 − k/N)
≤ 1 − 4δ(1 − δ).

*Step 4 (the two sides of (★★)).* Divide (★★) by N². By (ii) and (vi), the
left side is at most

  L := 1 − 4δ(1 − δ) + (√m + 2φ(τ) + ε₁)².

By Lemma 6 and ‖Re_T‖² = Np(1 − p)X, the right side equals f(p, b, ζ, X/N),
where

  f(p, β, γ, x) := (1 − 2p)² + (β(1 − 2p) + 2γ)² + 4p(1 − p)x.

For p ≤ 1/2, β ≥ 0, γ ≥ 0: ∂f/∂β = 2(1 − 2p)(β(1 − 2p) + 2γ) ≥ 0,
∂f/∂γ = 4(β(1 − 2p) + 2γ) ≥ 0 and ∂f/∂x = 4p(1 − p) ≥ 0. Moving one
coordinate at a time from (b₀, c₀, x₀) up to (b, ζ, X/N), along paths that stay
in {β ≥ 0, γ ≥ 0}, f does not decrease; so the right side is at least

  R := f(p, b₀, c₀, x₀).

At the limit point p = Q(θ), β = √m, γ = φ(θ), x = m, expanding the square
and using m(1 − 2p)² + 4p(1 − p)m = m and (1 − 2p)² = 1 − 4p(1 − p) gives

  f(Q(θ), √m, φ(θ), m) = 1 + m + 4[√m φ(θ)(1 − 2Q(θ)) + φ(θ)² − Q(θ)(1 − Q(θ))],

and, with L* := 1 − 4δ(1−δ) + (√m + 2φ(τ))²,

  f(Q(θ), √m, φ(θ), m) − L* = 4(Num − √m D) ≥ S.

*Step 5 (slack).* Error in L: L − L* = 2ε₁(√m + 2φ(τ)) + ε₁² ≤ 2ε₁(√a + 1) + ε₁².

Error in R. On the box K_a := [0, 1/2] × [0, √a] × [0, 1] × [−1, a + 1], the
partial derivatives of f are bounded by
|∂_p f| ≤ 4 + 4√a(√a + 2) + 4(a + 1), |∂_β f| ≤ 2(√a + 2),
|∂_γ f| ≤ 4(√a + 2), |∂_x f| ≤ 1, all at most K := 8 + 8a + 8√a. K_a is convex,
so f is K-Lipschitz on it for the ℓ¹ distance. Both (p, b₀, c₀, x₀) and
(Q(θ), √m, φ(θ), m) lie in K_a once 2√(aη) + 1/N ≤ 1 (note c₀ ≤ φ(θ) < 1 and
p ∈ [0, 1/2] by (6.1)). Their ℓ¹ distance is at most

  ε + √2 (aη)^{1/4} + 2ε + 1/N + 2√(aη),

using √m − b₀ ≤ √(2√(mη)) (from √u − √v ≤ √(u − v) for u ≥ v ≥ 0, and, when
m < 2√(mη), √m < √(2√(mη))). Therefore

  R − L ≥ S − K(3ε + √2(aη)^{1/4} + 2√(aη) + 1/N) − 2ε₁(√a + 1) − ε₁²,

uniformly in m ∈ (0, a]. Choose ε satisfying (6.1) and η ∈ (0, 1) so small
that

  K(3ε + √2(aη)^{1/4} + 2√(aη)) + 4ε(√a + 1) + 4ε² < S/2,

then N₀ ≥ max(N₁, N₂) so large that for N ≥ N₀: 2τ/N ≤ ε (so ε₁ ≤ 2ε and
2ε₁(√a + 1) + ε₁² ≤ 4ε(√a + 1) + 4ε²), K/N ≤ S/2, and 2√(aη) + 1/N ≤ 1.
Then R − L > S − S/2 − S/2 = 0, so (★★) holds on E, and with (i)
statement 1 holds with ρ = σ₁. ∎

Remarks.

- The probability bound holds for each N ≥ N₀, not only in the limit.
- θ < ∞ enters only through Q(θ) < δ (σ₁ inside the ball), the
  concentration of p, ζ, ‖Re_T‖, and the value of the limit; γ > 0 (σ₁ below ξ
  in energy) is *not* assumed — it holds at the optimum but the proof does not
  need it.
- T is defined from w = z + ν**1** rather than from z itself, because w is
  i.i.d. Gaussian, which makes p an i.i.d. Bernoulli average and ζ an
  i.i.d. nonnegative sum plus a scalar Gaussian term. Any
  other choice of T as a measurable function of (s, z, ν) independent of R
  would work in Lemmas 6 and 7; only Steps 2 and 3 use the specific form.

## 7. What the certification checks

[check_refstate.py](check_refstate.py) (standard library only) brackets π by Machin's formula,
1/√(2π) by rational square-root brackets, e^{−x²/2} and
∫_0^x e^{−u²/2}du by their Taylor series with geometric tail bounds (the
ratio of consecutive terms is non-increasing in absolute value, which
justifies the bound |remainder| ≤ |t_n|/(1 − |r_n|)), at x = τ₀ and x = θ₀.
It then uses monotonicity (x ↦ x(1 − x) increasing on (0, 1/2), φ(θ) ≥ 0,
1 − 2Q(θ) > 0) to bound Num and D from both sides and certifies
0 < τ₀ < θ₀, Q(θ₀) < Q(τ₀), 1 − 2Q(θ₀) > 0, Num > 0, D > 0 and
A₂(τ₀, θ₀) ≥ (Num_lo/D_hi)² ≥ 0.134716. An arb (python-flint) cross-check
of the same quantities agrees to all printed digits; it is not part of the
decision.

## 8. Gaps and limits (honest list)

- **No gap is known in the proof of Theorem 3**, but it has been checked only
  by its author (no independent review in this session).
- Lemma 4 (Sudakov–Fernique bound on E F) and Lemma 3 are imported from
  PROOF.md, not re-proved here.
- Citation details (theorem numbers in Laurent–Massart, Boucheron–Lugosi–
  Massart) are from memory and were not re-checked.
- Theorem 3 is about model G only. The ±1 model (model R) remains open for the
  same reason as in PROOF.md, Section 6 (hypothesis (U)); in addition, the
  decomposition of Lemma 5 uses rotation invariance of Gaussian rows and has
  no direct ±1 analogue.
- The supremum sup A₂ ≈ 0.1347166 is a floating-point value; only the value
  at the rational point (τ₀, θ₀) is certified.
- 0.134716 is a lower bound produced by this method; it is not claimed to be
  the capacity of model G.
