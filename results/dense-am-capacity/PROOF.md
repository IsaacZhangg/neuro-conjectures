# Proof: α_c(λ) = α_c^lb(λ) for Gaussian patterns, all λ ∈ (0, ∞]

This file proves the claim in [CLAIM.md](CLAIM.md) under the static retrieval
criterion stated there. Every assumption is listed in Section 0. Constants are
certified by `check.py` with interval arithmetic (python-flint / Arb).

## 0. Setting, notation, assumptions

- ξ^1, …, ξ^P iid N(0, I_N), P = ⌊e^{αN}⌋, α > 0 fixed, N → ∞.
- s_μ = ‖ξ^μ‖²/N. Its rate function (Cramér, χ²_N / N) is
  I(s) = (s − 1 − log s)/2 for s > 0. We use the exact Chernoff bounds
  P(s_μ ≤ t) ≤ e^{−N I(t)} for t ≤ 1 and P(s_μ ≥ t) ≤ e^{−N I(t)} for t ≥ 1,
  and the LDP lower bound liminf (1/N) log P(s_μ ∈ (a, b)) ≥ −inf_{(a,b)} I.
- Noise term, for λ ∈ (0, ∞): F^μ = (1/λ) log Σ_{ν≠μ} exp(λ ξ^μ·ξ^ν);
  F^μ_∞ = max_{ν≠μ} ξ^μ·ξ^ν. Pattern μ is retrieved iff N s_μ > F^μ.
  R_N = {all μ retrieved}; α_c(λ) = sup{α : P(R_N) → 1}.
- Facts used:
  (F1) Given ξ^μ, the overlaps ξ^μ·ξ^ν, ν ≠ μ, are iid N(0, N s_μ).
  (F2) ξ^μ = ‖ξ^μ‖ u^μ with u^μ uniform on the sphere and independent of
  ‖ξ^μ‖. For independent μ ≠ μ', ρ = u^μ·u^{μ'} has density
  c_N (1 − ρ²)^{(N−3)/2} on (−1, 1), c_N = Γ(N/2)/(Γ(1/2)Γ((N−1)/2)) ≤ √N,
  and ρ is independent of (s_μ, s_{μ'}).
  (F3) Monotonicity: adding patterns can only increase every F^μ (a new term
  in a log-sum-exp, or in a max), and adds new conditions. Coupling P ≤ P'
  on the same sequence, R_N(P') ⊂ R_N(P). So P(R_N) is nonincreasing in α at
  fixed N, and it suffices to prove the upper bound for α in a small interval
  just above the threshold.
  (F4) Mills: x φ(x)/(1 + x²) ≤ Φ̄(x) ≤ e^{−x²/2} for x > 0.

### The typical noise level and the threshold

For s > 0 let

    φ_typ(s) = α/λ + λ s/2   if λ² s < 2α,
    φ_typ(s) = √(2 α s)      if λ² s ≥ 2α,

(φ_typ(s) = √(2αs) for λ = ∞). This is the REM free energy: with
ψ_s(y) = (1/λ)(α − y²/(2s))_+ + y, one has sup_{y ≤ Y} ψ_s(y) ≤ max(φ_typ(s), Y)
for every Y (equality once Y is past the maximiser; only ≤ is used) (the concave part peaks at y = λ s when λ s < √(2αs); beyond
y = √(2αs), ψ_s(y) = y).

φ_typ(s) − s is positive at s → 0, and φ_typ is concave with slope < 1 once
it is below s, so there is a unique s0 = s0(α, λ) with φ_typ(s) > s for s < s0
and φ_typ(s) < s for s > s0. Explicitly s0 = 2α if λ ≥ 1 (and λ = ∞), and
s0 = α / c(λ) with c(λ) = λ(1 − λ/2) if λ < 1 (one checks λ² s0 < 2α ⇔ λ < 1).

Let J(t) = I(t) for t < 1 and J(t) = 0 for t ≥ 1, and define

    α_a(λ) = sup{ α : J(s0(α, λ)) > α },    α_b = (log 2)/4.

**Lemma 0 (closed form of the threshold).**
(i) α_c^lb(λ) = min(α_a(λ), α_b), where α_c^lb is the union-bound threshold of
LM eq. (10).
(ii) α_a(λ) = 1/(2e) for λ ≥ 1 (and λ = ∞), and for λ < 1,
α_a(λ) = −c W0(−a/e)/a with a = (1 − λ)², c = λ(1 − λ/2), W0 the principal
Lambert branch.
(iii) α_a is nondecreasing and continuous in λ. With
λ_x = inf{λ > 0 : α_a(λ) ≥ α_b}, α_c^lb(λ) = α_b for λ ≥ λ_x and
= α_a(λ) < α_b for λ < λ_x, and λ_x ∈ (0.700, 0.701) (certified).

*Proof.* LM's exponent is A(α, λ) = inf_s [I(s) + inf_{φ > s} I_{α,s,λ}(φ)]
(their eq. 9 with r² = s), where I_{α,s,λ}(φ) is +∞ for φ < φ_typ(s), 0 at
φ = φ_typ(s), and φ²/(2s) − α for φ > φ_typ(s) (their eq. 64). For s < s0 we
have φ_typ(s) > s, so the inner infimum is 0, attained at φ = φ_typ(s). For
s ≥ s0 we have φ_typ(s) ≤ s, so every admissible φ lies above φ_typ(s) and the
inner infimum is s/2 − α, approached as φ ↓ s (this is the single-overlap
branch, since φ > s ≥ φ_typ(s) ≥ √(2αs)). Hence

    A − α = min( J(s0) − α ,  inf_{s ≥ s0} [I(s) + s/2] − 2α ).

The map s ↦ I(s) + s/2 is convex with minimum log 2 / 2 at s = 1/2. If
s0 ≤ 1/2 the second term is log2/2 − 2α, positive iff α < α_b. If s0 > 1/2
it equals I(s0) + s0/2 − 2α ≥ I(s0) − α (use s0 ≥ 2α, which holds because
φ_typ(s0) = s0 and φ_typ(s) ≥ √(2αs)), so it never binds before the first
term. Moreover s0 > 1/2 and J(s0) > α force α < J(s0) < I(1/2) = 0.0966… < α_b
(I is decreasing on (0, 1)), so no α in (α_b, α_a) can enter the union-bound
set through this case, and {α : A − α > 0} is exactly the interval
(0, min(α_a, α_b)). This gives (i). For (ii): λ ≥ 1 gives I(2α) > α ⇔ −1 − log 2α > 0 ⇔
α < 1/(2e). For λ < 1 put x = α/c; I(x) = c x ⇔ a x − 1 = log x ⇔
x = −W(−a/e)/a, and the root below 1 is the W0 one. For (iii): φ_typ(s) is
nonincreasing in λ (d/dλ (α/λ + λs/2) = s/2 − α/λ² < 0 on λ²s < 2α), so s0 is
nonincreasing in λ, J(s0) nondecreasing, and the set {α : J(s0) > α} grows
with λ. So α_a is nondecreasing. It is continuous on (0, ∞] by (ii): W0 is
continuous on (−1/e, 0), and −c W0(−a/e)/a → c/e = 1/(2e) as λ ↑ 1 (a → 0,
W0(z) ~ z). With λ_x defined as the infimum in (iii), α_a < α_b below λ_x by
definition, and monotonicity with continuity gives α_a ≥ α_b from λ_x on.
`check.py` certifies α_a(0.700) < α_b < α_a(0.701); by continuity the first
inequality holds on a neighbourhood of 0.700 and the second on a
neighbourhood of 0.701, so 0.700 < λ_x < 0.701. ∎

The main theorem is then:

**Theorem.** For every λ ∈ (0, ∞]: if α < α_c^lb(λ) then P(R_N) → 1, and if
α > α_c^lb(λ) then P(R_N) → 0. Hence α_c(λ) = α_c^lb(λ); in particular
α_c(λ) = (log 2)/4 for all λ ≥ λ_x, including λ = ∞.

Proof: Section 1 (lower bound), Sections 2 and 3 (upper bound).

## 1. Lower bound: α < min(α_a, α_b) ⇒ P(R_N) → 1

By the union bound and exchangeability, P(R_N^c) ≤ P · P(F^1 ≥ N s_1). It
suffices to show P(F^1 ≥ N s_1) ≤ e^{−N(α + κ)} for some κ > 0.

Write X_ν = ξ^1·ξ^ν / N, ν = 2..P; given ξ^1 these are iid N(0, s/N) with
s = s_1 (F1), and let K = max_ν X_ν. The constants are chosen in this order.
First fix a compact s-window [σ, S] with 0 < σ < 1 < S such that
I(σ), I(S) > 2α + 1; then P(s_1 ∉ [σ, S]) ≤ 2 e^{−N(2α+1)}. Next take
ε, δ > 0 small enough that the error η = η(ε, δ) of Lemma 1 satisfies
η < σ/2, together with the further smallness conditions below. The cell
width h of the last step is chosen last.

**Lemma 1 (REM upper tail).** There are η(ε, δ) → 0 as ε, δ → 0 and an event
G with P(G^c | ξ^1) ≤ exp(−e^{εN/2}) for all N ≥ N0, uniformly in s ∈ [σ, S],
such that on G ∩ {K ≤ Y1}, F^1/N ≤ max(φ_typ(s), K) + η.

*Proof.* For λ = ∞ there is nothing to prove (F/N = K). For λ < ∞ take
levels y_k = −Y0 + kδ, k = 0, …, k_max, covering [−Y0, Y1], with
Y0 = (α + 2)/λ and Y1 = S + 1 (so Y1 > S ≥ s, which is all Section 1 needs).
Let n_k = #{ν : X_ν ∈ [y_k, y_{k+1})}, p_k = P(X ∈ [y_k, y_{k+1})) ≤
exp(−N d_k²/(2s)) with d_k the distance from 0 to the cell (F4), and
m_k = α − d_k²/(2s). Put L_k = ⌈exp(N(m_k^+ + ε))⌉. Since E n_k ≤ e^{N m_k},
the binomial tail P(n ≥ L) ≤ (e E n / L)^L for an integer L gives
P(n_k ≥ L_k) ≤ (e^{1 − εN})^{L_k} ≤ exp(−(εN − 1) e^{εN}) once εN ≥ 1. Let
G = {n_k < L_k for all k}, so that n_k ≤ exp(N(m_k^+ + ε)) on G. Its
complement is a union over O(1/δ) levels, so P(G^c | ξ^1) ≤ exp(−e^{εN/2})
for N ≥ N0, uniformly in s ∈ [σ, S].

Split Σ_ν e^{λ N X_ν} into the terms with X_ν < −Y0 (at most P e^{−λ N Y0} ≤
e^{−2N} in total), the terms with X_ν > Y1 (present only if K > Y1), and the
rest. On G ∩ {K ≤ Y1} the rest is at most

    M = (k_max + 1) · max_{k : y_k ≤ K} exp(N(m_k^+ + ε) + λ N y_{k+1})

(the maximum is empty only when K < −Y0, and then M = 0). So the whole sum is
at most e^{−2N} + M ≤ 2 max(e^{−2N}, M), and

    F/N ≤ (log 2)/(λN) + max( −2/λ , (1/(λN)) log M ).

For each level in the maximum there is a point y ≤ K, in the cell or at
distance at most δ from it, with

    (m_k^+ + ε)/λ + y_{k+1} ≤ ψ_s(y) + ε/λ + C δ.

Indeed, if y_k ≥ 0 the left side is ψ_s(y_k) + ε/λ + δ, and y = y_k; if
y_{k+1} ≤ 0 it is ψ_s(y_{k+1}) + ε/λ, and y = min(y_{k+1}, K); if the cell
contains 0 it is at most ψ_s(0) + ε/λ + δ, and y = min(0, K). Here C − 1 bounds the slope of ψ_s
on [−Y0, Y1 + δ], uniformly in s ∈ [σ, S] (|ψ_s'| ≤ (Y0 + Y1 + 1)/(λσ) + 1).
Hence

    (1/(λN)) log M ≤ sup_{y ≤ K} ψ_s(y) + ε/λ + C δ + log(k_max + 1)/(λN)
                   ≤ max(φ_typ(s), K) + ε/λ + C δ + log(k_max + 1)/(λN).

Since φ_typ(s) > 0 > −2/λ, the term −2/λ never exceeds max(φ_typ(s), K),
whatever the sign of K (this is what absorbs the e^{−2N} tail when K < 0).
Therefore F/N ≤ max(φ_typ(s), K) + η with η = ε/λ + (C + 1) δ, the extra δ
absorbing the two O(1/N) terms for N ≥ N0. η tends to 0 as ε, δ → 0,
uniformly for s ∈ [σ, S]. If K > Y1 we do not use the bound: that case is
inside the event {K ≥ s − η} treated below. ∎

Now on G ∩ {s_1 ∈ [σ, S]}, failure F^1 ≥ N s_1 forces (using Y1 > s_1 when
K > Y1) either φ_typ(s_1) ≥ s_1 − η or K ≥ s_1 − η. So

    P(F^1 ≥ N s_1) ≤ e^{−e^{εN/2}} + 2e^{−N(2α+1)} + P(φ_typ(s_1) ≥ s_1 − η)
                     + P(s_1 ∈ [σ, S], K ≥ s_1 − η).

*Short-pattern term.* φ_typ(s) ≥ s − η ⇔ s ≤ s0^η, where s0^η → s0 as η → 0
(φ_typ − s is strictly decreasing near s0). Since α < α_a, J(s0) > α; by
continuity pick η with J(s0^η) > α + 2κ. If s0^η < 1 the Chernoff bound gives
P(s_1 ≤ s0^η) ≤ e^{−N I(s0^η)} ≤ e^{−N(α + 2κ)}. (J(s0) > α > 0 forces s0 < 1.)

*Single-overlap term.* Given s_1 = s ∈ [σ, S] we have s − η > σ/2 > 0, so
P(K ≥ s − η) ≤ P Φ̄((s − η)√(N/s)) ≤ exp(N(α − (s − η)²/(2s))). Cover [σ, S]
by cells of width h; on each
cell P(s_1 ∈ cell) ≤ exp(−N inf_cell I) (Chernoff on the side away from 1; the
cell containing 1 has I ≥ 0). So the term is at most
(S/h) exp(N sup_{s ∈ [σ,S]} [α − I(s) − (s − η)²/(2s)] + N ω(h)) with ω(h) → 0.
As η, h → 0, inf_s [I(s) + (s − η)²/(2s)] → inf_s [I(s) + s/2] = (log 2)/2, so
the exponent tends to α − (log 2)/2 < −α − 2κ' because α < α_b.

Both terms are ≤ e^{−N(α+κ)} for small ε, δ, η, h, which proves the lower
bound. ∎

## 2. Upper bound, mechanism (b): α > α_b ⇒ P(R_N) → 0, every λ ∈ (0, ∞]

By (F3) it suffices to take α ∈ (α_b, α_*) with α_* = I(1/2) − h(1/2) − 1/4,
where

    h(r) = sup_{ρ ∈ (−1,1)} [ ½ log(1 − ρ²) − r/(1 + ρ) ].

`check.py` certifies h(1/2) ∈ −0.431446467943417… ± 10⁻⁵⁰ (maximizer
ρ* = [−(1 + r) + √((1 + r)² + 4r)]/2 ≈ 0.2808) and α_* = 0.2780200…, so
α_b = 0.17328… < α_*.

**Key observation.** If some ν ≠ μ has ξ^μ·ξ^ν > ‖ξ^μ‖² then μ fails, for
every λ: F^μ ≥ (1/λ) log e^{λ ξ^μ·ξ^ν} = ξ^μ·ξ^ν.

Fix r = 1/2 and small δ > 0. Split indices into A = {1, …, ⌊P/2⌋} and B (the
rest, |B| ≥ P/2). Let S_A = {μ ∈ A : s_μ ∈ (r − δ, r]}. Given the patterns in
A, the events E_ν = {∃ μ ∈ S_A : ξ^ν·ξ^μ > N s_μ}, ν ∈ B, are iid with
probability w(A). On R_N none of them occurs, so

    P(R_N) ≤ E[(1 − w(A))^{|B|}] ≤ E[exp(−|B| w(A))].

By Bonferroni, w(A) ≥ T1 − T2 with

    T1 = Σ_{μ ∈ S_A} Φ̄(√(N s_μ)),
    T2 = Σ_{μ < μ' ∈ S_A} P(Z_μ > √(N s_μ), Z_{μ'} > √(N s_{μ'}) | A),

where Z_μ = ξ·u^μ for a fresh ξ, so (Z_μ, Z_{μ'}) is a standard Gaussian pair
with correlation ρ_{μμ'} = u^μ·u^{μ'}.

*T1 is large.* q = P(s ∈ (r − δ, r]) satisfies e^{−N(I(r)+θ)} ≤ q ≤ e^{−N I(r)}
for N ≥ N0(θ), any θ > 0 (LDP lower bound on the open cell, I decreasing on
(0,1); Chernoff upper bound). |S_A| ~ Bin(⌊P/2⌋, q) has mean
m ≥ e^{N(α − I(r) − 2θ)}, and α − I(1/2) > α_b − I(1/2) = 0.0767… > 0, so
m → ∞ exponentially and P(|S_A| < m/2) ≤ e^{−m/8}. On {|S_A| ≥ m/2}, by (F4),

    T1 ≥ (m/2) Φ̄(√(N r)) ≥ e^{N(α − I(r) − r/2 − 3θ)}.

*T2 is small.* For thresholds a, a' ≥ b = √(N(r − δ)),
P(Z > a, Z' > a') ≤ P(Z + Z' > 2b) = Φ̄(b √(2/(1 + ρ))) ≤ e^{−N(r−δ)/(1+ρ)}.
By (F2), ρ is independent of the norms, so

    E T2 ≤ (P/2)² q² E[e^{−N(r−δ)/(1+ρ)}]
         ≤ e^{N(2α − 2I(r))} √N ∫ (1 − ρ²)^{(N−3)/2} e^{−N(r−δ)/(1+ρ)} dρ
         ≤ e^{N(2α − 2I(r) + h(r − δ) + θ)}

for N ≥ N0 (the integrand is ≤ exp(N[½ (1 − 3/N) log(1 − ρ²) − (r−δ)/(1+ρ)]),
whose sup converges to h(r − δ); h is continuous). By Markov,

    P(T2 > T1/2, |S_A| ≥ m/2) ≤ 2 e^{N(2α − 2I(r) + h(r−δ) + θ)} / e^{N(α − I(r) − r/2 − 3θ)}
                              = 2 e^{N(α − I(r) + h(r−δ) + r/2 + 4θ)}.

At r = 1/2 the exponent is α − α_* + (h(1/2 − δ) − h(1/2)) + 4θ < 0 for
δ, θ small, since α < α_*.

*Conclusion.* Outside an event of vanishing probability, w(A) ≥ T1/2 ≥
½ e^{N(α − I(1/2) − 1/4 − 3θ)} = ½ e^{N(α − (log 2)/2 − 3θ)}, because
I(1/2) + 1/4 = (log 2)/2. Then |B| w(A) ≥ ¼ e^{N(2α − (log 2)/2 − 3θ)} → ∞
since α > α_b = (log 2)/4 and θ is small. Hence P(R_N) → 0. ∎

Note the witness: a short pattern with ‖ξ^μ‖² ≈ N/2 (s ≈ 1/2), beaten by one
competitor whose overlap with it exceeds N/2, that is, whose component along
the pattern's direction u^μ is at least ‖ξ^μ‖ ≈ √(N/2). This is the same
optimum s = 1/2 that the union bound picks, which is why the bound is tight.

## 3. Upper bound, mechanism (a): α > α_a(λ) ⇒ P(R_N) → 0, for λ < 1

(For λ ≥ 1, α_a = 1/(2e) > α_b and Section 2 already covers every α > α_c^lb.)
By (F3) take α ∈ (α_a, α_a + ε) with ε small; then J(s0) < α. Also s0 < 1 at α = α_a (since
J(s0(α_a)) = α_a > 0) and s0 is continuous in α, so for ε small s0 < 1 and
I(s0) < α, and also λ² s0 < 2α (λ < 1). Pick t < s0 with I(t) < α and
λ² t < 2α (continuity), and small δ > 0.

Split A, B as in Section 2. The number of μ ∈ A with s_μ ∈ (t − δ, t] is
binomial with mean ≥ e^{N(α − I(t) − θ)} → ∞, so with probability → 1 there is
such a μ; let μ* be the smallest one (a function of A only). Condition on A.
By (F1), the overlaps ξ^{μ*}·ξ^ν, ν ∈ B, are iid N(0, N s) with s = s_{μ*}.
Let y = λ s and n = #{ν ∈ B : ξ^{μ*}·ξ^ν ≥ N y}. Then
p = Φ̄(λ√(N s)) ≥ e^{−N λ² s/2 − O(log N)} and E n = |B| p ≥
e^{N(α − λ² s/2) − O(log N)}, with α − λ² s/2 ≥ α − λ² t/2 > 0. By Chernoff,
n ≥ E n / 2 with probability ≥ 1 − e^{−E n/8}. On that event

    F^{μ*} ≥ (1/λ) log(n e^{λ N y}) ≥ N(α/λ − λ s/2 + λ s) − O(log N)
          = N φ_typ(s) − O(log N),

and φ_typ(s) − s ≥ φ_typ(t) − t > 0 for s ≤ t < s0 (φ_typ − s is decreasing
on (0, s0]: its slope is λ/2 − 1 < 0 in this branch). So F^{μ*} > N s_{μ*} for
N large: μ* fails. Hence P(R_N) → 0. ∎

## 4. Putting it together

Sections 2 and 3 give α_c(λ) ≤ min(α_a(λ), α_b) and Section 1 gives
α_c(λ) ≥ min(α_a(λ), α_b); by Lemma 0 this common value is α_c^lb(λ). In
fact the transition is sharp: P(R_N) → 1 below and → 0 above. ∎

## 5. Scope and limits

- The result is for the **static** criterion ‖ξ^μ‖² > F_λ^μ (LM eq. 6, and
  eq. 86 at λ = ∞), which is the event LM's union bound and α_c^lb are about.
  It does not address the dynamical definition (convergence of the iteration
  of LM eq. 2 from a nearby start). Remark: at x = ξ^μ the static event is
  exactly "the softmax weight on μ dominates"; turning that into convergence
  of the dynamics would need control of a neighbourhood of ξ^μ, not done here.
- Both mechanisms work for the threshold only through first- and
  second-moment estimates; no replica or cavity assumption is used.
- The value of λ_x is certified to an interval, not given in closed form.
- No Lean formalization was done (no toolchain in this environment; Mathlib
  coverage of the needed probability is thin). The constants are certified by
  interval arithmetic instead.

### Remark: a margin in the criterion does not move the threshold

*Added in the 2026-11-27 review (see REVIEW.md). It is not part of the
Theorem, it is a sketch, and it has not been independently reviewed.*

LM's criterion compares signal and noise on the exponential scale, so one can
ask whether the Theorem depends on using the exact inequality N s_μ > F^μ at
finite N. It does not. For γ ∈ ℝ say that μ is γ-retrieved if
N s_μ > F^μ + γN, and let R_N^γ be the event that every pattern is
γ-retrieved. A positive γ is stricter than the Theorem's criterion and a
negative γ is laxer; P(R_N^γ) is nonincreasing in γ, and (F3) holds for
every γ.

**Claim.** Fix λ ∈ (0, ∞]. For every α < α_c^lb(λ) there is γ0 > 0 with
P(R_N^γ) → 1 for all γ ≤ γ0. For every α > α_c^lb(λ) there is γ0 > 0 with
P(R_N^{−γ}) → 0 for all γ ≤ γ0. So the all-pattern threshold under the
γ-criterion tends to α_c^lb(λ) as γ → 0 from either side.

*Sketch.* Each section has room, and its exponents are continuous in γ.

- Section 1, stricter criterion (γ > 0). On G ∩ {s_1 ∈ [σ, S]}, failure
  F^1 ≥ N(s_1 − γ) forces φ_typ(s_1) ≥ s_1 − η − γ or K ≥ s_1 − η − γ. The
  rest of Section 1 goes through with η + γ in place of η, since it only uses
  that this quantity is small and below σ/2.
- Section 2, laxer criterion (margin −γ, γ > 0). A competitor with
  ξ^μ·ξ^ν > N(s_μ + γ) makes μ fail. Use these events in E_ν. The Gaussian
  thresholds become (s_μ + γ)√(N/s_μ). They are still ≥ b, so the bound on
  E T2 is unchanged. For s_μ ∈ (r − δ, r] they are ≤ √N (√r + γ/√(r − δ)),
  so T1 ≥ (m/2) Φ̄(√N (√r + γ/√(r − δ))), whose exponent tends to
  α − I(r) − r/2 as γ, θ → 0. The two strict inequalities the proof needs
  at r = 1/2 (T2 below T1/2, which uses α < α_*, and |B| T1 → ∞, which uses
  α > α_b) therefore survive for small γ.
- Section 3, laxer criterion. μ* fails by N(φ_typ(t) − t) − O(log N), a
  positive multiple of N, so every γ < φ_typ(t) − t is covered.
