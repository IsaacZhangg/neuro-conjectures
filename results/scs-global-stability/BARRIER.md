# Barrier: diagonal Lyapunov weights cannot pass g = 1/√2

PROOF.md proves global stability for g < 1/√2 with the Lyapunov function
Σ log cosh x_i, and notes that a diagonal weight, V_D = Σ d_i log cosh x_i,
works whenever J − I is Lyapunov diagonally stable (LDS). Whether such weights
reach past 1/√2 was left open there. This file proves that they do not: for
every fixed g > 1/√2, with probability tending to 1 exponentially fast, J − I
is not LDS, so no diagonal weight works. Any proof of the conjecture for
1/√2 < g < 1 therefore has to use something other than this family.

In the random-matrix and ecology literature this was an explicitly open
question. Gibbs, Grilli, Rogers and Allesina (Phys. Rev. E 98, 022410, 2018,
arXiv:1708.08837) write: "Under which conditions they are also diagonally
stable is an important open problem", and list as one possibility that
"negative definiteness is almost surely a necessary condition for diagonal
stability". Theorem 3 proves that possibility for zero-mean iid Gaussian
matrices, for each fixed g past the threshold (almost surely as N → ∞ along
any sequence, by Borel–Cantelli and the exponential bound), and the elliptic
remark below extends it to zero-mean Gaussian correlated pairs (J_ij, J_ji).
Non-Gaussian entries and nonzero mean, which their question also covers, are
not treated. Our search (NOVELTY.md) found no later paper that settles it.

This is a theorem about the method, not about the network. It does not show
that the network fails to be globally stable for g > 1/√2.

## Setting and definitions

As in PROOF.md: G is N×N with iid N(0,1) entries (diagonal included),
J = (g/√N) G, A = J − I. For a matrix M, M_s = (M + Mᵀ)/2 and
M_a = (M − Mᵀ)/2. A real matrix A is *Lyapunov diagonally stable* if there is
a diagonal D with positive entries such that DA + AᵀD is negative definite.
This is exactly the condition under which the argument of PROOF.md, Theorem 1,
goes through with V_D in place of V.

## Theorem 3

Fix g > 1/√2. There are constants c > 0 and N₀, depending only on g, such
that for all N ≥ N₀

    P(J − I is Lyapunov diagonally stable) ≤ exp(−cN).

**Corollary.** Under the same conditions, with probability at least
1 − exp(−cN), for every positive diagonal D the function
V_D(x) = Σ_i d_i log cosh x_i is not a Lyapunov function for the flow: in
every neighbourhood of 0 there are points x ≠ 0 where dV_D/dt > 0. (Proved
after Theorem 3.)

Together with PROOF.md, Theorem 2, this pins the asymptotic threshold for
diagonal stability of (g/√N)G − I at exactly g = 1/√2. Equivalently,
inf_D λ_max of the symmetric part of D^{1/2}(G/√N)D^{−1/2}, over positive
diagonal D, tends to √2 in probability. At finite N the best diagonal weight
lowers this quantity a little below λ_max of the symmetric part of G/√N
(PROOF.md, "Diagonal weights", measured at g = 1), and the theorem says the
gain disappears as N → ∞.

## Proof

**Lemma 4 (the easy direction of Barker, Berman and Plemmons 1978).** If A is
Lyapunov diagonally stable, then for every positive semidefinite X ≠ 0 some
diagonal entry of AX is negative.

*Proof.* Let D be the diagonal weight and M = DA + AᵀD ≺ 0. Since X ⪰ 0 and
X ≠ 0, tr(MX) = tr(X^{1/2} M X^{1/2}) < 0. By cyclicity
tr(AᵀDX) = tr((AᵀDX)ᵀ) = tr(XDA) = tr(DAX), so
tr(MX) = 2 tr(DAX) = 2 Σ_i d_i (AX)_ii. If every (AX)_ii were ≥ 0 this would
be ≥ 0. ∎

So it is enough to build, with high probability, a positive semidefinite
X ≠ 0 with (AX)_ii ≥ 0 for every i.

**Two facts about Gaussian G.**

(F1) G_s and G_a are independent. Each pair (G_ij + G_ji, G_ij − G_ji) for
i < j is a pair of independent Gaussians, the diagonal of G goes only into
G_s, and different pairs are independent.

(F2) G_s is the Gaussian orthogonal ensemble: its diagonal entries are N(0,1),
its entries above the diagonal are N(0,1/2), all independent, so its density
is proportional to exp(−tr H²/2). Hence U G_s Uᵀ has the same law as G_s for
every fixed orthogonal U. The entries of G_a above the diagonal are iid
N(0, 1/2).

Write H = J_s = (g/√N) G_s and K = J_a = (g/√N) G_a, so that A = H − I + K,
H and K are independent, and K_ij for i < j are iid N(0, g²/(2N)).

**The certificate.** Fix a with 1 < a < √2 g and let δ = (√2 g − a)/2. Let
P = 1_{[a,∞)}(H) be the orthogonal projection onto the span of the
eigenvectors of H with eigenvalue ≥ a, and m = rank P. Take X = P. Then

    (AX)_ii = ((H − I)P)_ii + (KP)_ii.

Writing P = Σ_{k: λ_k ≥ a} v_k v_kᵀ with orthonormal eigenvectors v_k,
((H − I)P)_ii = Σ_k (λ_k − 1) v_{k,i}² ≥ (a − 1) P_ii. And
(KP)_ii = Σ_{j≠i} K_ij P_ji =: T_i (K has zero diagonal). So

    (AX)_ii ≥ (a − 1) P_ii + T_i,                                    (2)

and it suffices to show that, with probability ≥ 1 − exp(−cN), m ≥ 1 and
|T_i| ≤ (a − 1) P_ii / 2 for every i. (Half the margin is kept so that the
corollary gets a strict inequality.)

**Step 1: many eigenvalues above a.** Let h(x) = min(1, (x − a)₊/δ), so
0 ≤ h ≤ 1_{[a,∞)} and h is (1/δ)-Lipschitz, and let F = (1/N) Σ_k h(λ_k(H)),
so m ≥ NF. As a function of G with the Frobenius norm, F is
(g/(δN))-Lipschitz: by the Hoffman–Wielandt inequality,
Σ_k |λ_k(H) − λ_k(H′)| ≤ √N ‖H − H′‖_F ≤ √N (g/√N) ‖G − G′‖_F. Gaussian
concentration (as in PROOF.md, Theorem 2(b)) gives

    P(F ≤ E F − t) ≤ exp(−δ² N² t² / (2 g²)).

By Wigner's semicircle law the mean empirical spectral measure of H converges
weakly to the semicircle law σ of radius √2 g. Let 2η := ∫ h dσ, which is
> 0 because h is positive on (a, √2 g]. So E F ≥ 3η/2 for N ≥ N₀, and with
t = η/2,

    P(m < ηN) ≤ exp(−δ² η² N² / (8 g²)).                             (3)

**Step 2: the diagonal of P.** Let U be a Haar-distributed orthogonal matrix
independent of G. By (F2), H and UHUᵀ have the same law, and
1_{[a,∞)}(UHUᵀ) = U P Uᵀ with the same rank m. So (m, P) and (m, UPUᵀ) have
the same law. Conditional on (m, P), UPUᵀ is the projection onto a uniformly
random m-dimensional subspace, so conditional on m, P is such a projection.
For such a projection P_ii = |Pe_i|² has the law of
Q = (y_1² + … + y_m²)/(y_1² + … + y_N²) with y iid N(0,1). By the χ² tail
bounds of Laurent and Massart (2000), with x = m/16 the numerator is at least
m − 2√(mx) = m/2 and the denominator at most N + 2√(Nx) + 2x ≤ 2N, except
with probability at most 2 exp(−m/16). So, given m ≥ 1,

    P(P_ii < m/(4N) | m) ≤ 2 exp(−m/16)   for each i.                (4)

**Step 3: the antisymmetric noise.** K is independent of H, hence of P. Given
P, T_i = Σ_{j≠i} K_ij P_ji is a sum of independent centred Gaussians, so it is
N(0, σ_i²) with

    σ_i² = (g²/(2N)) Σ_{j≠i} P_ij² ≤ (g²/(2N)) P_ii,

because P² = P gives Σ_j P_ij² = P_ii. The Gaussian tail gives

    P(|T_i| > (a − 1) P_ii / 2 | P) ≤ 2 exp(−(a − 1)² P_ii² / (8σ_i²))
                                   ≤ 2 exp(−(a − 1)² N P_ii / (4g²)).

On the event P_ii ≥ m/(4N) this is at most 2 exp(−(a − 1)² m / (16g²)).

**Step 4: union bound.** On the event m ≥ ηN, by (4) and Step 3 and a union
bound over i, the probability that some i has P_ii < m/(4N) or
|T_i| > (a − 1) P_ii / 2 is at most

    2N exp(−ηN/16) + 2N exp(−(a − 1)² η N / (16g²)).

Adding (3), the bad event has probability at most exp(−cN) for N ≥ N₀ with,
for example, c = (1/2) min(η/16, (a − 1)² η/(16g²)). Off the bad event, m ≥ 1,
P ⪰ 0, P ≠ 0 and by (2)

    (AP)_ii ≥ (a − 1) P_ii / 2 ≥ (a − 1) m / (8N) > 0   for every i,     (5)

so by Lemma 4 A = J − I is not Lyapunov diagonally stable. ∎

**Proof of the corollary.** Off the bad event, let D be any positive diagonal
matrix and M_D = D(J − I) + (J − I)ᵀD. As in Lemma 4,
tr(M_D P) = 2 Σ_i d_i (AP)_ii > 0 by (5), so M_D has an eigenvalue μ > 0
with a unit eigenvector u. Along the flow,

    dV_D/dt = φᵀ D(−x + Jφ) = (1/2) φᵀ M_D φ + Σ_i d_i φ_i (φ_i − x_i),

with φ = tanh x. At x = ru with r > 0: |φ − ru| ≤ r³/3 (since
|tanh s − s| ≤ |s|³/3) and |φ_i| ≤ r|u_i|, so (1/2) φᵀM_Dφ ≥ (r²/2) μ − C r⁴
and Σ d_i φ_i(φ_i − x_i) ≥ −max_i d_i r⁴/3, with C depending on D and J.
Hence dV_D/dt ≥ (r²/2) μ − C′r⁴ > 0 at x = ru for all small r > 0. ∎

## Remarks

- **Exponential rate, constants not optimised.** The proof gives some c > 0;
  it does not try to make c or N₀ explicit or good.
- **Without self-couplings.** If G has zero diagonal, the symmetric part is
  H′ = H − diag(H). Replace H by the GOE matrix H″ = H′ + diag(ξ) with ξ an
  independent N(0, g²/N) vector, and take X = 1_{[a,∞)}(H″). Then (2) becomes
  (AX)_ii ≥ (a − 1 − |ξ_i|) P_ii + T_i. P depends on ξ, but only the bound
  max_i |ξ_i| ≤ (a − 1)/4 is used, and it fails with probability at most
  2N exp(−(a − 1)² N/(32 g²)). Running Steps 3 and 4 with margin
  (a − 1)/4 in place of (a − 1)/2 gives the same conclusion with probability
  at least 1 − exp(−c′N).
- **Elliptic couplings.** Suppose instead that the pairs (G_ij, G_ji), i < j,
  are iid centred Gaussian pairs with unit variances and correlation
  ρ ∈ (−1, 1), independent of the diagonal (any distribution with
  max_i |G_ii| ≤ √(4 log N) w.h.p.). Then G_ij + G_ji and G_ij − G_ji are
  still independent, the off-diagonal part of H has iid N(0, g²(1 + ρ)/(2N))
  entries, and K has iid N(0, g²(1 − ρ)/(2N)) entries above the diagonal.
  Repeating the no-self-coupling argument (an auxiliary independent diagonal
  ξ ~ N(0, g²(1 + ρ)/N) makes the matrix GOE; the margin must now absorb
  |J_ii| + |ξ_i|, both O(√(log N / N)) w.h.p.), J − I fails to be LDS with
  probability tending to 1 for every fixed g > 1/√(2(1 + ρ)), exponentially
  fast when the diagonal is Gaussian. The negative-definiteness
  threshold is g < 1/√(2(1 + ρ)), the constant in Clenet, El Ferchichi and
  Najim (arXiv:2205.15591, Prop. 2.3), so for every ρ ∈ (−1, 1) the two
  thresholds coincide. (ρ = 1 is the symmetric case, where they coincide
  trivially.) Nonzero mean interactions are not covered.
- **Other activations.** Lemma 4 and Theorem 3 are about the matrix J − I
  only, so the barrier applies to the Matsuoka / Forti–Tesi criterion for any
  slope-restricted activation, and to every argument that reduces to Lyapunov
  diagonal stability of J − I.
- **Random Lotka–Volterra.** Clenet, Massol and Najim (2023) get a unique,
  globally stable equilibrium for random Lotka–Volterra systems from the
  Takeuchi–Adachi criterion, that is, from Lyapunov diagonal stability of
  I − Γ/(α√N), with the same √2 threshold. For Gaussian
  Γ with iid N(0,1) entries, Theorem 3 with g = 1/α says that this route stops
  exactly at α = √2. We did not check their non-Gaussian settings.
- **What the theorem does not say.** It says nothing about whether the origin
  is globally stable for 1/√2 < g < 1. Lyapunov functions that are not of the
  diagonal form, or arguments that use the dynamics, are not ruled out.

## Numerical check (method fixed before running)

barrier_check.py builds, for N in {500, 1000, 2000, 4000} and
g in {0.8, 0.9, 1.0}, one Gaussian matrix each (seed 2026), and two
candidate certificates:

1. X = P = 1_{[a,∞)}(J_s) with a = (1 + √2 g)/2, as in the proof;
2. X = Σ_k (λ_k − 1)₊ v_k v_kᵀ, which uses every eigenvalue above 1.

It records m (the rank), min_i (AX)_ii / mean_i (AX)_ii, the fraction of
negative diagonal entries, and whether X certifies non-LDS (all entries ≥ 0).
It records the same with the constants of the proof: min_i P_ii · N/m and
max_i |T_i| / ((a − 1) P_ii). The theorem is asymptotic, so at small N the
certificate is expected to fail; the check only looks for the fraction of
negative entries falling as N grows. Output goes to barrier-check-output.json.
This check is a guide, not part of the proof.

## References

- G. P. Barker, A. Berman, R. J. Plemmons, "Positive diagonal solutions to the
  Lyapunov equations", *Linear and Multilinear Algebra* 5 (1978) 249–256.
- B. Laurent, P. Massart, "Adaptive estimation of a quadratic functional by
  model selection", *Ann. Statist.* 28 (2000) 1302–1338, Lemma 1.
- G. W. Anderson, A. Guionnet, O. Zeitouni, *An Introduction to Random
  Matrices*, CUP 2010 (semicircle law; Hoffman–Wielandt, Lemma 2.1.19).
- S. Boucheron, G. Lugosi, P. Massart, *Concentration Inequalities*, OUP 2013,
  Thm 5.6.
- T. Gibbs, J. Grilli, T. Rogers, S. Allesina, "Effect of population
  abundances on the stability of large random ecosystems", *Phys. Rev. E* 98
  (2018) 022410, arXiv:1708.08837, Section VII.
- M. Clenet, H. El Ferchichi, J. Najim, arXiv:2205.15591, Prop. 2.3.
- M. Clenet, F. Massol, J. Najim, *J. Math. Biol.* (2023), arXiv:2205.00735,
  Theorem 2 (via the Takeuchi–Adachi criterion).
