# Proof

This file proves Theorem 1 (Gaussian interference patterns) completely, and
reduces the standard ±1 model to one unproved hypothesis (U) in Theorem 2.
Section 6 explains why (U) is the real obstacle. Nothing here is claimed to be
new mathematics: the constant is the first-level value of Stojnic (2024), and
the method is the Sudakov–Fernique inequality. What this file adds is a
self-contained statement and proof with every step written out, a certified
constant, and an exact statement of what is still missing for ±1 patterns.

## 1. Setting and notation

N ≥ 2 neurons, spins σ ∈ {−1, +1}^N. The target pattern is ξ ∈ {−1, +1}^N.
There are M ≥ 0 further ("interference") patterns, the rows g_1, …, g_M of an
M × N matrix G, so P = M + 1 patterns in total and the load is α = P/N.

Hebbian couplings, no self-coupling:

  J_ij = (1/N) (ξ_i ξ_j + Σ_μ G_μi G_μj) for i ≠ j,  J_ii = 0.

Local field h_i(σ) = Σ_{j≠i} J_ij σ_j. Energy H(σ) = −(1/2) Σ_{i≠j} J_ij σ_i σ_j.
Hamming distance d(σ, τ) = #{i : σ_i ≠ τ_i}. A configuration σ is a **stable
fixed point** if σ_i h_i(σ) > 0 for every i (no single spin flip, and no step
of asynchronous zero-temperature dynamics σ_i ← sign h_i(σ), changes σ).

- **Model G (Theorem 1):** the entries of G are independent N(0, 1).
- **Model R (Theorem 2, the standard Hopfield model):** every pattern,
  target included, has independent uniform ±1 entries, and the target is any
  one of the P patterns.

Φ and φ are the standard normal distribution function and density,
Q(t) = 1 − Φ(t). For τ > 0 put

  δ(τ) = Q(τ),  A(τ) = ((δ(1−δ) − φ(τ)²) / φ(τ))²  if δ(1−δ) > φ(τ)², else 0.

Σ_k is the set of σ ∈ {±1}^N with exactly k coordinates equal to −1, and
**1** is the all-ones vector. ‖·‖ is the Euclidean norm.

## 2. Two elementary lemmas

**Lemma 1 (energy identity).** For every σ,

  H(σ) = −(1/(2N)) [ (ξ·σ)² + ‖Gσ‖² ] + (1/(2N)) [ N + Σ_{μ,i} G_μi² ].

*Proof.* For any vector a, Σ_{i≠j} a_i a_j σ_i σ_j = (a·σ)² − Σ_i a_i², because
σ_i² = 1. Apply this to a = ξ and to a = g_μ for each μ and sum. ∎

**Lemma 2 (barrier implies a fixed point).** Fix k ≥ 1. Suppose that
H(σ) > H(ξ) for every σ with d(σ, ξ) = k, and that H(σ) ≠ H(σ′) whenever
σ′ is obtained from σ by flipping one spin. Then there is a stable fixed
point σ* with d(σ*, ξ) ≤ k − 1. In particular ξ·σ*/N ≥ 1 − 2(k−1)/N.

*Proof.* The ball B = {σ : d(σ, ξ) ≤ k − 1} is finite and contains ξ. Let σ*
minimize H on B. Let σ′ be σ* with spin i flipped. Then d(σ′, ξ) ≤ k. If
d(σ′, ξ) ≤ k − 1, then H(σ′) ≥ H(σ*) by minimality, and H(σ′) ≠ H(σ*) by
assumption, so H(σ′) > H(σ*). If d(σ′, ξ) = k, then H(σ′) > H(ξ) ≥ H(σ*).
Flipping spin i changes the energy by exactly H(σ′) − H(σ*) = 2 σ*_i h_i(σ*)
(the terms of H containing i sum to −σ_i h_i(σ), and the flip reverses their
sign). So σ*_i h_i(σ*) > 0 for every i. ∎

**Lemma 2′ (ties).** Without the no-ties assumption, the same proof gives a
σ* in B with σ*_i h_i(σ*) ≥ 0 for every i: a fixed point of the dynamics
under the convention that a spin with zero field keeps its value.

**Lemma 3 (no ties in model G).** In model G with N ≥ 2, almost surely no two
configurations that differ in exactly one spin have equal energy.

*Proof.* Let σ′ differ from σ in one spin. By Lemma 1,
2N (H(σ′) − H(σ)) = (ξ·σ)² − (ξ·σ′)² + Σ_μ g_μᵀ (σσᵀ − σ′σ′ᵀ) g_μ, a polynomial
in the entries of G. Its quadratic part is zero only if σσᵀ = σ′σ′ᵀ, that is
σ′ = ±σ, which is impossible for a single flip when N ≥ 2. A nonzero
polynomial vanishes only on a Lebesgue-null set, and there are finitely many
pairs. ∎

## 3. Reduction to a Gaussian maximum

**Gauge.** In model G, replacing σ_i by ξ_i σ_i and G_μi by G_μi ξ_i leaves the
law of G unchanged and maps ξ to **1**, preserving energies and Hamming
distances. So assume ξ = **1**.

For σ ∈ Σ_k we have **1**·σ = N − 2k, so by Lemma 1

  H(σ) > H(**1**)  ⟺  ‖Gσ‖² < ‖G**1**‖² + N² − (N − 2k)² = ‖G**1**‖² + 4k(N − k).

So the barrier at radius k holds as soon as

  (★)  F² < ‖G**1**‖² + 4k(N − k),  where F := max_{σ∈Σ_k} ‖Gσ‖.

## 4. The expected maximum

**Theorem A (Sudakov–Fernique).** Let (X_t) and (Y_t), t in a finite set T,
be centered jointly Gaussian processes with E(X_s − X_t)² ≤ E(Y_s − Y_t)² for
all s, t. Then E max_T X ≤ E max_T Y. (Fernique 1975; Vitale 2000;
Chatterjee 2005 gives a short proof; Adler and Taylor, *Random Fields and
Geometry*, Theorem 2.2.3.)

**Lemma 4.** E F ≤ √(NM) + 2Nφ(t) + 2t(k − NQ(t)) for every t ∈ ℝ. In
particular, if k = ⌈δ(τ)N⌉ then E F ≤ √(NM) + 2Nφ(τ) + 2τ.

*Proof.* Index set T = Σ_k × S^{M−1}. Put X_{σ,y} = yᵀGσ and
Y_{σ,y} = √N g·y + h·σ with g ∈ ℝ^M, h ∈ ℝ^N independent standard Gaussian
vectors. Using ‖σ‖² = N and ‖y‖ = 1,

  E(X_{σ,y} − X_{σ′,y′})² = ‖yσᵀ − y′σ′ᵀ‖²_F = 2N − 2(y·y′)(σ·σ′),
  E(Y_{σ,y} − Y_{σ′,y′})² = N‖y − y′‖² + ‖σ − σ′‖² = 4N − 2N(y·y′) − 2(σ·σ′),

and the difference is 2(N − σ·σ′)(1 − y·y′) ≥ 0. For any finite T′ ⊂ T,
Theorem A gives E max_{T′} X ≤ E max_{T′} Y ≤ E sup_T Y. Take finite T′_n
increasing to Σ_k × D with D a countable dense subset of the sphere. Since
y ↦ yᵀGσ is continuous, max_{T′_n} X increases to sup_T X = F, and
max_{T′_1} X is integrable, so by monotone convergence E F ≤ E sup_T Y.

Now sup_T Y = √N ‖g‖ + max_{σ∈Σ_k} h·σ and E‖g‖ ≤ (E‖g‖²)^{1/2} = √M. For
σ ∈ Σ_k with minus set S, h·σ = Σ_i h_i + 2 Σ_{i∈S} (−h_i). For any reals
a_1, …, a_N, any t, and any |S| = k, Σ_{i∈S} a_i ≤ kt + Σ_i (a_i − t)_+.
With a_i = −h_i (standard normal) and E(Z − t)_+ = φ(t) − tQ(t),

  E max_{σ∈Σ_k} h·σ ≤ 0 + 2[kt + N(φ(t) − tQ(t))] = 2Nφ(t) + 2t(k − NQ(t)).

With t = τ and k = ⌈Q(τ)N⌉ we have 0 ≤ k − NQ(τ) < 1. ∎

## 5. Theorem 1 and its proof

**Theorem 1.** In model G, fix τ > 0 with A(τ) > 0, put δ = δ(τ) and
k = ⌈δN⌉, and let M = M(N) satisfy limsup M/N < A(τ). There are c > 0 and
N₀ such that for N ≥ N₀, with probability at least 1 − 2e^{−cN},

1. H(σ) > H(ξ) for every σ with d(σ, ξ) = k, and
2. there is a stable fixed point σ* with d(σ*, ξ) ≤ k − 1, hence
   ξ·σ*/N ≥ 1 − 2δ − 2/N.

**Corollary.** At τ₀ = 1.8259 (δ ≈ 0.0339), A(τ₀) ≥ 0.12948 (certified by
check.py). So for any load α < 0.12948 the target pattern of model G has an
energy barrier and a stable fixed point with overlap at least 0.932 with
probability 1 − e^{−cN}. Numerically, sup_τ A(τ) = 0.1294900 at
τ ≈ 1.82587; this supremum is not needed for the Corollary and is not
certified.

*Proof of Theorem 1.* By the gauge, take ξ = **1**. Statement 2 follows from
statement 1 by Lemmas 2 and 3 (on the probability-one event of Lemma 3).
Choose a with limsup M/N < a < A(τ) and N₁ with M ≤ aN for N ≥ N₁.

Since 4φ√m + 4φ² is increasing in m and equals 4δ(1−δ) at m = A(τ) (by the
definition of A), the slack s := 4δ(1−δ) − 4φ√a − 4φ² is positive.

*Concentration of F.* For matrices G, G′,
|F(G) − F(G′)| ≤ max_{σ∈Σ_k} ‖(G − G′)σ‖ ≤ ‖G − G′‖_op √N ≤ ‖G − G′‖_F √N,
so F is √N-Lipschitz in the Frobenius norm. By Gaussian concentration for
Lipschitz functions (Tsirelson–Ibragimov–Sudakov; e.g. Boucheron, Lugosi and
Massart, *Concentration Inequalities*, Theorem 5.6),
P(F ≥ EF + εN) ≤ exp(−ε²N/2).

*Lower tail of ‖G**1**‖².* G**1** has independent N(0, N) entries, so
‖G**1**‖²/N is χ² with M degrees of freedom. By Laurent and Massart (2000,
Lemma 1), P(χ²_M ≤ M − 2√(Mx)) ≤ e^{−x}; take x = ηN.

*Combination.* Off an event of probability at most e^{−ε²N/2} + e^{−ηN},
Lemma 4 gives, with m = M/N ≤ a,

  F/N ≤ √m + 2φ + ε + 2τ/N  and  ‖G**1**‖²/N² ≥ m − 2√(mη).

Also k/N ∈ [δ, δ + 1/N] ⊂ [δ, 1 − δ] for large N, and x ↦ 4x(1−x) is at least
4δ(1−δ) on [δ, 1 − δ], so 4k(N−k)/N² ≥ 4δ(1−δ). Hence (★) holds if

  (√m + 2φ + ε′)² < m − 2√(mη) + 4δ(1−δ),  with ε′ = ε + 2τ/N,

that is, if 4φ√m + 4φ² + 2ε′(√m + 2φ) + ε′² + 2√(mη) < 4δ(1−δ). The left side
is at most 4φ√a + 4φ² + 2ε′(√a + 2φ) + ε′² + 2√(aη). Choose ε, η > 0 so small,
and N₀ ≥ N₁ so large, that 2ε′(√a + 2φ) + ε′² + 2√(aη) < s for N ≥ N₀. Then
(★) holds, which is statement 1. Take c = min(ε²/2, η). ∎

Remarks.

- The probability bound holds for each N ≥ N₀, not only in the limit.
- Theorem 1 is about one target pattern with Gaussian interference. In
  model G the interference patterns are not ±1 configurations, so "every
  stored pattern" has no meaning there.
- The supremum A* = 0.1294900 agrees with Stojnic's first-level value
  α_c(NLT,1) ≈ 0.129490. His reported closed form
  erf(x)²/(2x²) − 1 + erf(x)², evaluated at x = τ*/√2 ≈ 1.291085, matches
  A(τ*) to about 3·10⁻⁹ in floating point (check.py prints both). Whether the
  two expressions agree exactly at the optimum was not checked symbolically.
  The agreement is expected: his first level of lifting is this comparison.

## 6. Binary patterns: Theorem 2 and the missing step

In model R, gauge by the target: σ_i ↦ ξ_i σ_i and X_μi := ξ^μ_i ξ_i for the
other M = P − 1 patterns. X has independent uniform ±1 entries independent of
ξ, and Lemma 1 holds with G replaced by X (now Σ X_μi² = MN is a constant).
The reduction (★) is unchanged.

**Hypothesis (U).** For k = ⌈δ(τ)N⌉ and M/N → α,

  E max_{σ∈Σ_k} ‖Xσ‖ ≤ √(NM) + 2Nφ(τ) + o(N)  for X with ±1 entries.

**Theorem 2 (conditional).** Assume (U) for some τ with α < A(τ). Then in
model R, with probability tending to one (exponentially fast), every one of
the P stored patterns has an energy barrier at radius ⌈δN⌉ and, by Lemma 2′,
a fixed point within distance ⌈δN⌉ − 1 (with the zero-field convention).

*Proof sketch, all steps standard.* F(X) = max_{σ∈Σ_k} ‖Xσ‖ is convex in X
and √N-Lipschitz, and the entries are independent and bounded, so
Talagrand's concentration inequality for convex Lipschitz functions gives
P(|F − med F| ≥ r) ≤ 4 exp(−r²/(16N)) (Ledoux, *The Concentration of Measure
Phenomenon*, Corollary 4.10, rescaled from [0, 1] to [−1, 1]); the mean and
median differ by O(√N). ‖X**1**‖² = Σ_μ S_μ² with S_μ = Σ_i X_μi independent,
E S_μ² = N and E S_μ⁴ ≤ 3N², so the Chernoff bound with e^{−λy} ≤ 1 − λy + λ²y²/2
gives P(‖X**1**‖² ≤ (1 − ε)MN) ≤ exp(−ε²M/6). The combination step of
Theorem 1 then applies verbatim with (U) in place of Lemma 4, giving failure
probability e^{−cN} per pattern, and a union bound over the P = O(N)
patterns. ∎

**Why (U) is not proved here.** The Sudakov–Fernique inequality needs
Gaussian processes. The standard routes from ±1 to Gaussian entries do not
reach (U):

- *Contraction.* Writing G = X ∘ |G| and using convexity gives
  E max ‖Xσ‖ ≤ √(π/2) E max ‖Gσ‖. The factor 1.2533 destroys the bound.
- *Entry-wise Lindeberg replacement with a soft maximum.* The error terms
  are of order β^j Σ_{μ,i} ⟨|u_μ|^{j+1}⟩ where u = Xσ/‖Xσ‖ under the Gibbs
  measure, so they are small only if u is delocalized. On Σ_k a single row
  can have |x_μ·σ| up to about 2kN/N ≈ 2δN, the same order as ‖Xσ‖ ≈ 0.5N,
  so delocalization would have to be proved for near-maximizers, uniformly
  along the interpolation. Choosing a fixed temperature does not help: the
  entropy cost H(δ)/β and the Lindeberg error are both of order one per spin.
- *Recent rigorous work points the same way.* Barra and Durante (2026,
  arXiv:2610.41938, unrefereed, proofs written with AI assistance) prove that
  replacing all binary patterns by Gaussian ones gives only a one-sided
  bound, and that the Rademacher ground-state gain dominates the Gaussian
  one (their Corollary 9.3). Their collective-replacement error,
  (13/3)β²ϱδ²e^{4√β δ} per spin, grows with β and is not usable at zero
  temperature at this precision. Ali, Boursier, Chen, Ho, Issa and Mourrat
  (2026, arXiv:2611.03306, unrefereed) show the Bernoulli free energy equals
  a Gaussian one with finitely many extra Bernoulli fields, again because the
  law of the patterns matters. The direction that would be needed for (U),
  binary at most Gaussian, is exactly the one these papers say can fail
  without extra structure.

What would close the gap: a proof that on Σ_k (overlap at least 1 − 2δ with
the target) no other pattern can carry a macroscopic share of ‖Xσ‖ at a
near-maximizer, strong enough to run a Lindeberg or collective-replacement
argument at zero temperature with error o(N). On Σ_k every non-target overlap
is at most 2δ + O(√(log N / N)) ≈ 0.068, which is the structure such a proof
would use.

## 7. References for cited theorems

- X. Fernique (1975), Régularité des trajectoires des fonctions aléatoires
  gaussiennes. LNM 480. R. A. Vitale (2000), Some comparisons for Gaussian
  processes, Proc. AMS 128, 3043–3046. S. Chatterjee (2005), An error bound in
  the Sudakov–Fernique inequality, arXiv:math/0510424.
- B. Laurent and P. Massart (2000), Adaptive estimation of a quadratic
  functional by model selection, Ann. Statist. 28, Lemma 1.
- S. Boucheron, G. Lugosi, P. Massart (2013), Concentration Inequalities,
  Theorem 5.6 (Gaussian concentration).
- M. Ledoux (2001), The Concentration of Measure Phenomenon, Corollary 4.10.

Page and theorem numbers in this section were taken from memory of the
standard texts and were not re-checked against the books in this session.
