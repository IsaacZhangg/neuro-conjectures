# Proof: the K = 3 tree committee machine stores every pattern up to α = 1.55

Agent-drafted (Claude). Not yet independently reviewed.

## 0. Model and statement

**Tree committee machine (TCM), K = 3.** There are n = 3N inputs split into three
disjoint blocks of N. Weights are w = (w₁, w₂, w₃) with w_k ∈ ℝᴺ. A pattern is
ξ = (ξ₁, ξ₂, ξ₃), ξ_k ∈ ℝᴺ, with output

  σ(ξ) = sign( Σ_k sign⟨w_k, ξ_k⟩ ).

We are given m = ⌊α n⌋ patterns ξ^μ with all entries i.i.d. N(0,1) and labels
y_μ ∈ {±1} that are fixed, or random and independent of the patterns. The
patterns are *stored* by w if σ(ξ^μ) = y_μ for every μ. The (storage) capacity
α_c(3) is the threshold in α of the event "some w stores all m patterns" as
N → ∞.

Known before this note (see SOURCES.md): the best proved upper bound is
α_c(3) ≤ 3.43 (Stojnic 2023, partially lifted random duality theory); replica
calculations give about 3.0 to 3.14 (one step of replica symmetry breaking) and
simulations 2 to 2.75. We found no explicit proved lower bound in the literature.

> **Theorem 1.** For every α ≤ 31/20 = 1.55,
> P( some w stores all ⌊αn⌋ patterns ) → 1 as N → ∞. Hence α_c(3) ≥ 1.55, in the
> sense lim inf_{N→∞} P(storable at load α) = 1 for every α ≤ 1.55.

The constant comes from the three-stage construction below. Every numerical
inequality it needs is checked in Arb ball arithmetic by [check.py](check.py)
(section 6).

**Remark (elementary bound 4/3).** Using only Cover's theorem one gets α_c(3) ≥ 4/3:
branch 1 is made positive on 2N(1−δ) chosen patterns; the rest are positive on
branch 1 with probability 1/2 each, independently, because those rows are
independent of w₁. Branch 2 is made positive on the branch-1 failures plus
further chosen patterns up to 2N(1−δ), and branch 3 on all patterns that failed
once. The load count gives 3α < 4. Theorem 1 improves this with convex
optimisation in stages 1 and 2.

## 1. Reduction

Put x^μ_k = y_μ ξ^μ_k. These are i.i.d. N(0, I_N) (labels independent of
patterns and ±1, Gaussian symmetric). Since sign is odd, σ(ξ^μ) = y_μ iff
Σ_k sign⟨w_k, x^μ_k⟩ > 0. It therefore suffices to find w such that for every
μ **at least two** of the three numbers ⟨w_k, x^μ_k⟩ are strictly positive (then
the sum of signs is at least 2 − 1 = 1). Write X_k ∈ ℝ^{m×N} for the matrix with
rows x^μ_k; X₁, X₂, X₃ are independent. Let β = m/N → 3α (we use β = 3α = 4.65 at
α = 1.55; a slightly smaller load only helps, see the end of section 5).

Fixed constants (exact rationals, as in check.py):

  κ = 1/1000, ρ₁ = 59/250, ρ₂ = 7/25, a = 975253/10⁶.

## 2. The construction

For ρ ∈ (0,1) define the convex losses

  ℓ_R(h) = (1−ρ)(κ − h)₊ − ρh  ("soft"),  ℓ_H(h) = ι_{[κ,∞)}(h) − ρh  ("hard"),

where ι is the convex indicator (0 on the set, +∞ off it).

**Stage 1.** Let ŵ₁ be any minimiser of Σ_{μ=1}^m ℓ_R(⟨x^μ₁, w⟩) over ‖w‖ ≤ 1, with
ρ = ρ₁. Let E₁ = {μ : ⟨x^μ₁, ŵ₁⟩ ≤ 0}.

**Stage 2.** On the event |E₁| ≤ aN, let H ⊇ E₁ be E₁ together with the smallest
indices outside E₁ so that |H| = ⌊aN⌋, and R = [m] ∖ H. (H, R are functions of X₁
only.) Let ŵ₂ be any minimiser over ‖w‖ ≤ 1 of
Σ_{μ∈R} ℓ_R(⟨x^μ₂, w⟩) + Σ_{μ∈H} ℓ_H(⟨x^μ₂, w⟩) with ρ = ρ₂, that is: minimise
(1−ρ₂)Σ_R(κ − h_μ)₊ − ρ₂ Σ_{all μ} h_μ subject to h_μ ≥ κ on H. Let
E₂ = {μ ∈ R : ⟨x^μ₂, ŵ₂⟩ ≤ 0}.

**Stage 3.** Let T = H ∪ E₂ (a function of X₁, X₂ only) and let ŵ₃ be any vector
with ⟨x^μ₃, ŵ₃⟩ > 0 for all μ ∈ T, if one exists.

**Why this stores everything.** Take any μ. If μ ∉ E₁, branch 1 is positive.
If μ ∈ H, branch 2 is ≥ κ > 0. If μ ∈ R ∖ E₂, branch 2 is positive. If μ ∈ T,
branch 3 is positive. A pattern in E₁ lies in H ⊆ T (branches 2, 3); a pattern
in H ∖ E₁ has branches 1, 2; a pattern in R ∖ E₂ has branches 1, 2; a pattern in
E₂ ⊆ R has branches 1, 3. So every pattern has two positive branches.

It remains to show that, with probability → 1: stage 2 is feasible and well
defined, and |T| ≤ (2 − δ)N for some δ > 0 so that stage 3 succeeds.

## 3. Lemma 1 (Cover, Wendel)

> Let x₁, …, x_p be i.i.d. N(0, I_N) and p ≤ (2 − δ)N. Then
> P(∃ w: ⟨w, x_i⟩ > 0 ∀ i) = 2^{1−p} Σ_{j=0}^{N−1} C(p−1, j) ≥ 1 − e^{−c(δ)N}.

The identity is Wendel's theorem (1962) for symmetric distributions in general
position (Cover 1965 counts the same dichotomies). The sum is
P(Bin(p−1, 1/2) ≤ N−1), and N − 1 ≥ (p−1)/2 + (δ/4)N for large N, so Hoeffding
gives the bound. Because T depends only on (X₁, X₂), which is independent of X₃,
Lemma 1 applies conditionally on T.

## 4. Lemma 2 (law of the margins at a convex optimum)

**Set-up.** Let H, R be disjoint deterministic index sets with |H|/N → L_H ∈ [0,1)
and |R|/N → L_R > 0, X ∈ ℝ^{(|H|+|R|)×N} with i.i.d. N(0,1) entries, κ > 0,
ρ ∈ (0,1), and consider

  (P)  min_{‖w‖≤1} F(w) := Σ_{μ∈R} ℓ_R(h_μ) − ρ Σ_{μ∈H} h_μ  s.t. h_μ ≥ κ (μ ∈ H),  h = Xw.

For x > 0 and G ~ N(0,1) define the *displacements*

  D^R_x(g) = max(xρ, min(x, κ − g)),  D^H_x(g) = max(xρ, κ − g),

and

  c(x) = L_H E[D^H_x(G)²] + L_R E[D^R_x(G)²],
  P(x) = L_H E[(G + D^H_x)D^H_x] + L_R E[(G + D^R_x)D^R_x].

Note g + D^R_x(g) = prox_{xℓ_R}(g) and g + D^H_x(g) = prox_{xℓ_H}(g): indeed
−∂ℓ_R(h) = {1−ρ+ρ} = {1} for h < κ, [ρ, 1] at κ, {ρ} for h > κ, and the prox
condition is prox(g) − g ∈ −x ∂ℓ(prox(g)); similarly for ℓ_H, whose
subdifferential at κ is (−∞, −ρ]. Each D_x(g) is non-decreasing in x for every
g, so c is continuous and strictly increasing, with c(0+) = L_H E[(κ−G)₊²] and
c(∞) = ∞.

> **Lemma 2.** Assume (i) c(x*) = 1 for some x* > 0 (unique by monotonicity),
> (ii) P(x*) > 0, and (iii) if L_H > 0: 2κ √L_H < 1 − √L_H. Then for every ε > 0,
> with probability → 1, problem (P) is feasible and **every** minimiser ŵ satisfies
>
>   #{μ ∈ R : ⟨x_μ, ŵ⟩ ≤ 0} ≤ (Φ(κ − x*) + ε)|R|.

Stage 1 is Lemma 2 with H = ∅, L_R = β. Stage 2 is Lemma 2 conditionally on X₁,
with L_H = a, L_R = β − a. Condition (iii) at L_H = a = 0.975253 reads
0.00198 < 0.01245, true.

### Proof of Lemma 2

Write m' = |H| + |R|, Lip for Lipschitz constants, and work on the event 𝓔_N
(probability → 1, standard Gaussian matrix bounds, Davidson and Szarek 2001):
‖X‖_op ≤ (1 + √β + δ)√N; if L_H > 0, the smallest singular value of the
|H| × N block X_H is ≥ (1 − √L_H − δ)√N; ‖z‖² / N ∈ [1−δ, 1+δ] and
max_μ |g_μ| ≤ 2√(log N) for the Gaussian vectors z ∈ ℝᴺ, g ∈ ℝ^{m'} used below.
Here δ > 0 is small and fixed.

**Step 1: feasibility and a Slater point.** If L_H > 0, let
w₀ = X_Hᵀ(X_H X_Hᵀ)⁻¹(2κ𝟙). Then X_H w₀ = 2κ𝟙 and
‖w₀‖ ≤ 2κ√|H| / s_min(X_H) ≤ 2κ√L_H/(1 − √L_H − δ) + o(1) < 1 by (iii).
So (P) is feasible, with every hard constraint slack by κ at w₀.

**Step 2: bounded multipliers and an exact penalty.** (P) is a finite-dimensional
convex program satisfying Slater's condition, so a Lagrange multiplier
u* ∈ ℝ^H_{≥0} for the constraints κ − h_μ ≤ 0 exists, and
F* = min_{‖w‖≤1} [F(w) + Σ_H u*_μ(κ − ⟨x_μ, w⟩)] ≤ F(w₀) − κ‖u*‖₁.
F is a sum of 1-Lipschitz functions of the coordinates of h, so
F(w₀) − F* ≤ ‖X(w₀ − ŵ)‖₁ ≤ √m' ‖X‖_op · 2 ≤ C₀N on 𝓔_N. Hence
‖u*‖_∞ ≤ ‖u*‖₁ ≤ C₀N/κ. Fix M_N = C₀N/κ + 1 and the Lipschitz loss
ℓ_H^M(h) = M_N(κ − h)₊ − ρh. By the exact penalty theorem for convex programs
(e.g. Bertsekas, *Convex Optimization Theory*, Prop. 5.5.1: if M exceeds
‖u*‖_∞, the ℓ₁-penalised problem has the same minimisers), the minimisers of (P)
are exactly the minimisers of

  (P_M)  min_{‖w‖≤1} Σ_R ℓ_R(h_μ) + Σ_H ℓ_H^M(h_μ),  h = Xw.

All losses in (P_M) are convex and L_N := M_N + 1 Lipschitz.

**Step 3: lifting and the Gaussian min-max theorem.** Let
𝓛(h) = Σ_R ℓ_R(h_μ) + Σ_H ℓ_H^M(h_μ), C₁ = 10(1 + √β), K_N = L_N √m' + 1, and
B_h = {‖h‖ ≤ C₁√N}. For any closed S ⊆ ℝ^{m'} define

  Φ_N(S) = min_{‖w‖≤1, h ∈ S∩B_h} max_{‖u‖≤K_N} uᵀXw − uᵀh + 𝓛(h)
      = min_{‖w‖≤1, h∈S∩B_h} 𝓛(h) + K_N‖Xw − h‖.

Because 𝓛 is L_N√m'-Lipschitz in ℓ₂ and K_N exceeds that, on 𝓔_N (where
Xw ∈ B_h) Φ_N(ℝ^{m'}) = min_{‖w‖≤1} 𝓛(Xw) = value of (P_M), and for every S,
Φ_N(S) ≤ min{𝓛(Xw) : ‖w‖≤1, Xw ∈ S} (take h = Xw).

The convex Gaussian min-max theorem (Thrampoulidis, Oymak and Hassibi 2015,
Theorem 3; built on Gordon 1988) compares Φ_N with

  φ_N(S) = min_{‖w‖≤1, h∈S∩B_h} max_{‖u‖≤K_N} ‖w‖gᵀu + ‖u‖zᵀw − uᵀh + 𝓛(h),

with g ∈ ℝ^{m'}, z ∈ ℝᴺ independent standard Gaussian vectors: for every t,
(a) P(Φ_N(S) < t) ≤ 2P(φ_N(S) ≤ t) for any closed S, and (b)
P(Φ_N(ℝ^{m'}) > t) ≤ 2P(φ_N(ℝ^{m'}) ≥ t), because the sets are convex and compact
and the objective is convex-concave. Maximising over the direction of u and then
minimising over the direction of w at fixed r = ‖w‖ (z ᵀw ≥ −r‖z‖):

  φ_N(S) = min_{r∈[0,1], h∈S∩B_h} 𝓛(h) + K_N(‖h − rg‖ − r‖z‖)₊.

**Step 4: the auxiliary problem.** Let V_N(r) = min{𝓛(h) : ‖h − rg‖ ≤ r‖z‖}.
Projecting h onto that ball changes 𝓛 by at most L_N√m' times the distance,
so φ_N(ℝ^{m'}) = min_r V_N(r) on 𝓔_N (the minimisers lie inside B_h, since
‖h‖ ≤ r‖g‖ + r‖z‖ ≤ (√β + 1 + 2δ)√N). By Lagrange duality for this
separable convex problem with one ball constraint,

  V_N(r)/N = max_{λ≥0} { (1/N) Σ_μ e_μ(r g_μ; 1/λ) − λ r²‖z‖²/(2N) },
  e_μ(v; x) = min_h ℓ_μ(h) + (h − v)²/(2x)  (Moreau envelope).

On 𝓔_N, |g_μ| ≤ 2√log N, so the prox of ℓ_H^M equals the prox of ℓ_H at every
coordinate (the penalty slope M_N is never reached), and the law of large
numbers gives, uniformly for r in [0,1] and λ in compacts of (0,∞),

  V_N(r)/N → v(r) := max_{λ≥0} d(r, λ),
  d(r, λ) = L_H E[e_H(rG; 1/λ)] + L_R E[e_R(rG; 1/λ)] − λr²/2.

Using ∂_v e(v; x) = −D_x(v)/x and ∂_λ e = D²/2 (with x = 1/λ, D = prox − v):
∂_λ d(1, λ) = (c(1/λ) − 1)/2, so at r = 1 the dual optimum is λ* = 1/x*, and

  ∂_r d(r, λ*)|_{r=1} = −(1/x*)(L_H E[G D^H] + L_R E[G D^R] + 1) = −P(x*)/x*  < 0

by (i) and (ii) (since 1 = c(x*)). Here we used x*-scaling: prox_{xℓ}(rg) for
the losses above, at r = 1, is g + D_x(g).

*v is convex*: the constraint set {(r, h): ‖h − rg‖ ≤ r‖z‖} is a convex cone, so
V_N is convex in r, and so is the limit v. Weak duality gives v(r) ≥ d(r, λ*)
for all r with equality at r = 1. Hence v(1−η₀) ≥ v(1) + P(x*)η₀/(2x*) for some
small η₀ > 0, and by convexity v(r) ≥ v(1) + c₂(1 − r) for all r ∈ [0,1], with
c₂ = P(x*)/(2x*) > 0. So the auxiliary minimum is attained at r = 1 in the
limit.

**Step 5: deviation.** Let ψ(h) = min(1, max(0, 1 − h/κ)), which is
(1/κ)-Lipschitz, equals 1 for h ≤ 0 and 0 for h ≥ κ. For the soft coordinates at
r = 1, ψ(G + D^R_{x*}(G)) = 0 unless G < κ − x* (for κ − x* ≤ G ≤ κ − x*ρ the prox
is exactly κ, above it is > κ). So e* := E ψ(G + D^R_{x*}(G)) ≤ Φ(κ − x*). Let

  S_ε = { h : (1/|R|) Σ_{μ∈R} ψ(h_μ) ≥ e* + ε }.

We show φ_N(S_ε) ≥ N(v(1) + γ) w.h.p. for some γ = γ(ε) > 0. Take (r, h) with
h ∈ S_ε ∩ B_h.

* If dist(h, ball_r) > η√N (ball_r = {‖h − rg‖ ≤ r‖z‖}), the penalty term is
  ≥ K_N η√N ≥ η N^{3/2}, far above N(v(1) + 1), because 𝓛 ≥ −ρ‖h‖₁ ≥ −C N on B_h.
* Otherwise let h' be the projection of h onto ball_r. Then 𝓛(h) + penalty ≥
  𝓛(h'), and since ψ is (1/κ)-Lipschitz, h' ∈ S_{ε − η'} with η' = η√N·√|R|/(κ|R|),
  which is small when η is small.
  + If r ≤ 1 − η₁: 𝓛(h') ≥ V_N(r) ≥ N(v(1) + c₂η₁ − o(1)).
  + If r ∈ [1 − η₁, 1]: for every feasible h', with x* = 1/λ*,
    𝓛(h') ≥ 𝓛(h') + (λ*/2)(‖h' − rg‖² − r²‖z‖²)
        ≥ N d_N(r, λ*) + (λ*/2)‖h' − h_r‖²,
    where h_r = (prox_{x*ℓ_μ}(r g_μ))_μ is the unique minimiser of the
    λ*-strongly convex function 𝓛(h) + (λ*/2)‖h − rg‖², and d_N → d uniformly.
    By the law of large numbers and since the prox is 1-Lipschitz in r,
    (1/|R|)Σ_R ψ(h_{r,μ}) ≤ e* + ε/4 for η₁ small and N large; with
    h' ∈ S_{ε−η'} and ψ (1/κ)-Lipschitz, ‖h' − h_r‖ ≥ κ(ε/2)√|R| (for η' ≤ ε/4).
    Since d(r, λ*) ≥ d(1, λ*) − ω(η₁) = v(1) − ω(η₁) with ω(η₁) → 0,
    𝓛(h') ≥ N(v(1) − ω(η₁) − o(1) + λ*κ²ε²L_R/8).

Choosing η, η₁ small gives γ > 0. By Step 3(a),
Φ_N(S_ε) ≥ N(v(1) + γ/2) w.h.p.; by Step 3(b) and Step 4,
Φ_N(ℝ^{m'}) ≤ N(v(1) + γ/4) w.h.p. If a minimiser ŵ of (P_M) had Xŵ ∈ S_ε, then
Φ_N(S_ε) ≤ 𝓛(Xŵ) = Φ_N(ℝ^{m'}), a contradiction. So every minimiser has
(1/|R|)Σ_R ψ(⟨x_μ, ŵ⟩) < e* + ε, and #{μ ∈ R : ⟨x_μ,ŵ⟩ ≤ 0} ≤ Σ_R ψ(⟨x_μ,ŵ⟩)
< (Φ(κ − x*) + ε)|R|. By Step 2 the same holds for (P). ∎

## 5. Proof of Theorem 1

Fix ε > 0 small.

*Stage 1.* Lemma 2 with L_H = 0, L_R = β = 4.65, ρ = ρ₁. check.py certifies a
bracket x₁* ∈ [x₁⁻, x₁⁺] = [0.808353…, 0.808354…] for the root of c, P > 0 on the
whole bracket (P ≥ 1.86·10⁻⁴), and β Φ(κ − x₁⁻) < a. Hence, w.h.p.,
|E₁| ≤ (Φ(κ − x₁*) + ε)m ≤ aN for ε small, so H and R are defined, with
|H|/N → a and |R|/N → β − a.

*Stage 2.* Conditionally on X₁, apply Lemma 2 with L_H = a, L_R = β − a,
ρ = ρ₂. check.py certifies the bracket [0.587626…, 0.587627…] for x₂*, P > 0
on it (P ≥ 1.27·10⁻³), and condition (iii) holds. Hence w.h.p. stage 2 is
feasible and |E₂| ≤ (Φ(κ − x₂*) + ε)|R|.

*Stage 3.* |T| ≤ |H| + |E₂| ≤ N[a + (β − a)(Φ(κ − x₂⁻) + ε)] + o(N). check.py
certifies a + (β − a)Φ(κ − x₂⁻) = 1.999505… < 2. Taking ε small, |T| ≤ (2 − δ)N for
some δ > 0, and Lemma 1 (conditionally on X₁, X₂) gives ŵ₃.

All three events have probability → 1, so the probability that w = (ŵ₁, ŵ₂, ŵ₃)
stores all patterns tends to 1 at α = 1.55. For α < 1.55, add dummy patterns up
to load 1.55 (storing more patterns stores the original ones). ∎

## 6. What check.py certifies

All quantities are closed-form in Φ and the Gaussian density φ:
E[D²] and E[GD] for the two displacement shapes reduce to
∫(c − g)²φ(g)dg over intervals, and E[G D] = −P(D is not at its cap) by
Gaussian integration by parts (D is absolutely continuous with D' = −1 on the
middle region, 0 elsewhere). check.py evaluates them in Arb ball arithmetic at
200 bits, finds the brackets for x* by bisection on the certified sign of
c(x) − 1, encloses P on the whole bracket with one ball evaluation, and checks
the six inequalities used above (including condition (iii) of Lemma 2). Output: [check-output.json](check-output.json).

## 7. Limits

* Gaussian patterns only. For ±1 inputs a universality argument would be
  needed; Cover's theorem is distribution-free but Lemma 2 is not.
* This is a lower bound from one explicit construction. The gap to the proved
  upper bound 3.43 stays wide; the construction is greedy (each stage optimises
  its own error count) and convex (stage optimisers cannot reach the non-convex
  optimum). Simulations in the literature put the true capacity at 2 to 2.75.
* The finite-size simulation [sim_pipeline.py](sim_pipeline.py) is a sanity
  check of the limit formulas only (it stored all patterns at α = 1.4, N = 250,
  and α = 1.5, N = 400, one seed each); it is not part of the proof.
