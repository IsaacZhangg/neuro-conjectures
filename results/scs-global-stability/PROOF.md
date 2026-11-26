# Proof: global stability of the SCS network for g < 1/√2

## Setting

    dx/dt = F(x) := -x + J φ(x),   φ(x) = (tanh x_1, ..., tanh x_N),

with J = (g/√N) G and G an N×N matrix of iid N(0,1) entries (diagonal
included). Write A_s = (A + Aᵀ)/2 for the symmetric part of a matrix and
λ_max for the largest eigenvalue of a symmetric matrix. Let

    W = G_s / √N = (G + Gᵀ) / (2√N),   so that   J_s = g W.

"Globally asymptotically stable" (GAS) means: the origin is Lyapunov stable and
every solution converges to 0.

## Theorem 1 (deterministic)

Let J be any real N×N matrix with λ := λ_max(J_s) < 1. Then the origin is GAS.
Moreover, with V(x) = Σ_i log cosh x_i, every solution satisfies

    V(x(t)) ≤ V(x(0)) · exp(-2(1-λ) κ t),   κ = (tanh M / M)²,   M = V(x(0)) + log 2,

so convergence is exponential with a rate that depends only on λ and on V(x(0)).

This is the case D = I of the classical criterion "J − I Lyapunov diagonally
stable ⇒ global stability" for x' = −x + Jφ(x) with slope-restricted φ
(Matsuoka 1992; Forti and Tesi 1995). We give the short proof to keep the
result self-contained.

**Proof.** F is globally Lipschitz (constant 1 + ‖J‖), so solutions exist and
are unique for all t ≥ 0.

*Properties of V.* For s ∈ ℝ: log cosh s ≥ 0 with equality only at s = 0;
log cosh s ≤ s²/2 (since (log cosh)′(s) = tanh s ≤ s on s ≥ 0, and evenness); and log cosh s ≥ |s| − log 2 (since cosh s ≥ e^{|s|}/2). So V is
C¹, positive definite and radially unbounded, with |x_i| ≤ V(x) + log 2 for
every i and 2V(x) ≤ |x|².

*Derivative.* ∇V(x) = φ(x), so along solutions

    dV/dt = φ·(−x + Jφ) = −φ·x + φᵀ J φ = −φ·x + φᵀ J_s φ.

For each i, tanh x_i has the sign of x_i and |tanh x_i| ≤ |x_i|, so
φ_i x_i ≥ φ_i². Hence φ·x ≥ |φ|² and

    dV/dt ≤ −|φ|² + λ|φ|² = −(1 − λ)|φ|²,                          (1)

which is < 0 whenever x ≠ 0 (φ = 0 only at x = 0). A C¹ function that is
positive definite and radially unbounded with negative definite derivative
proves GAS (Lyapunov's direct method with Barbashin–Krasovskii radial
unboundedness; e.g. Khalil, *Nonlinear Systems*, 3rd ed., Thm 4.2).

*Rate.* V is non-increasing, so |x_i(t)| ≤ V(x(0)) + log 2 = M for all i and
t ≥ 0. Since tanh is concave on [0, ∞), |tanh s| ≥ (tanh M / M)|s| for
|s| ≤ M, so |φ|² ≥ κ|x|² ≥ 2κV. With (1), dV/dt ≤ −2(1−λ)κ V, and Grönwall
gives the bound. ∎

## Theorem 2 (random coupling)

(a) E λ_max(W) ≤ √2 for every N.

(b) For every N and s ≥ 0, P(λ_max(W) ≥ √2 + s) ≤ exp(−N s²/2).

(c) Consequently, for every N and s > 0, with probability at least
1 − exp(−N s²/2), the origin is GAS **simultaneously for every
0 ≤ g ≤ 1/(√2 + s)** (same G, coupling J = (g/√N)G). In particular, for fixed
0 < g < 1/√2,

    P(origin is GAS) ≥ 1 − exp(−N (1 − √2 g)² / (2 g²)) → 1.

**Proof of (a).** λ_max(W) = sup_{u ∈ S^{N−1}} X_u with
X_u = uᵀ G u / √N = N^{−1/2} Σ_{i,j} G_ij u_i u_j, a centred Gaussian process
with continuous paths on the sphere. For unit vectors u, v with c = u·v,

    E(X_u − X_v)² = N^{−1} Σ_{i,j} (u_i u_j − v_i v_j)² = N^{−1} |uuᵀ − vvᵀ|_F² = (2 − 2c²)/N.

Let Y_u = √(2/N) ξ·u with ξ ~ N(0, I_N). Then E(Y_u − Y_v)² = (2/N)|u − v|² =
(4 − 4c)/N, and (4 − 4c) − (2 − 2c²) = 2(1 − c)² ≥ 0. By the Sudakov–Fernique
inequality (e.g. Adler and Taylor, *Random Fields and Geometry*, Thm 2.2.3;
apply it on finite subsets of the sphere and take an increasing union, using
path continuity),

    E sup_u X_u ≤ E sup_u Y_u = √(2/N) E|ξ| ≤ √(2/N) · √N = √2.

**Proof of (b).** As a function of G ∈ ℝ^{N×N} with the Frobenius norm,
G ↦ λ_max(W) is (1/√N)-Lipschitz: |λ_max(A) − λ_max(B)| ≤ ‖A − B‖_op ≤
‖A − B‖_F for symmetric A, B (Weyl), and ‖(H + Hᵀ)/2‖_F ≤ ‖H‖_F. Gaussian
concentration for Lipschitz functions (Tsirelson–Ibragimov–Sudakov; e.g.
Boucheron, Lugosi and Massart, *Concentration Inequalities*, Thm 5.6) gives
P(λ_max(W) ≥ E λ_max(W) + s) ≤ exp(−N s²/2), and (a) finishes it.

**Proof of (c).** On the event λ_max(W) < √2 + s we have
λ_max(J_s) = g λ_max(W) < g(√2 + s) ≤ 1 for every 0 < g ≤ 1/(√2 + s), and
Theorem 1 applies (g = 0 is trivial). For fixed 0 < g < 1/√2 take
s = 1/g − √2 > 0, so g = 1/(√2 + s) and s² = (1 − √2 g)²/g². (Only g ≥ 0
is meant throughout: g is the gain, and for g < 0 one would need the event
for λ_max(−W) instead.) ∎

## How far this method can go

- **The constant is sharp for the identity weight.** λ_max(W) → √2 almost
  surely (W is a Wigner matrix with off-diagonal variance 1/(2N);
  Bai and Yin 1988). So the Lyapunov function Σ log cosh x_i certifies
  exactly the range g < 1/√2 asymptotically: near 0,
  dV/dt ≈ −|x|² + xᵀJ_s x, which is positive along the top eigenvector of
  J_s once λ_max(J_s) > 1, so V fails to be a Lyapunov function there.
- **Diagonal weights.** V_D = Σ d_i log cosh x_i works whenever
  D(J − I) + (J − I)ᵀD ≺ 0 (Lyapunov diagonal stability of J − I), by the
  same computation. Numerically (check.py, check 4) the best diagonal D found
  at g = 1 lowers λ_max of sym(D^{1/2} J D^{−1/2}) by about 0.09 at N = 125,
  0.08 at N = 250, 0.06 at N = 500 and 0.055 at N = 1000 (mean of two
  matrices each), and the spread of the optimal log d shrinks with N
  (0.34, 0.28, 0.25, 0.22). We do not know whether the gain survives as N → ∞.
  The asymptotic threshold in g for Lyapunov diagonal stability of
  J − I = (g/√N)G − I is,
  as far as we found, an open random-matrix question; it lies in
  [1/√2, 1] since diagonal stability forces the spectral abscissa of gG
  below 1. The same open question governs the random Lotka–Volterra model
  (see SOURCES.md; we did not check that this
  equivalence is stated there in exactly this form).
- **Multiplier methods (not proved here, a remark).** Lur'e–Postnikov
  functions x ᵀPx + Σ d_i ∫φ with S-procedure constraints, and frequency-domain
  Zames–Falb multipliers for repeated odd monotone nonlinearities, at
  frequency 0 need a static multiplier M₀ in the doubly-dominant class with
  M₀(I − J) + (I − J)ᵀM₀ᵀ ≻ 0. For dense iid J with entries of size N^{−1/2},
  diagonal dominance limits the off-diagonal part of M₀ to row sums of the
  size of its diagonal, which does not change the edge of the spectrum at
  leading order (though a sparse off-diagonal part of size O(1) is not
  obviously negligible). We therefore expect these methods to stop at the diagonal
  threshold above. Going past it likely needs facts specific to tanh (for
  example tanh′(s) < 1 for s ≠ 0, and s − tanh s ≥ tanh³(s)/3) or arguments
  that use the dynamics. This remark is a heuristic, not a theorem.

## What remains open

Global asymptotic stability with high probability for 1/√2 ≤ g < 1. The
origin is locally exponentially stable there (spectral radius of J tends to
g, no outliers; Geman 1986 and Bai and Yin 1986 for the spectral radius,
Tao 2013 for outliers; Bordenave and Chafaï's survey), and simulations
from wide initial conditions converge to 0 (check.py, check 3), but no
global argument is proved.

## Assumptions, stated once

Real Gaussian couplings with variance g²/N, self-couplings included (they are
harmless: removing the diagonal changes W by a diagonal matrix of norm
about √(2 log N / N), and the proof of (a) goes through with Σ_{i≠j},
since then E(X_u − X_v)² = (2 − 2c² − Σ_i (u_i² − v_i²)²)/N ≤ (2 − 2c²)/N,
giving the same bound; (b) is unchanged). Activation tanh. The GAS part of Theorem 1 (not the
explicit rate) holds for any locally Lipschitz φ with sφ(s) ≥ φ(s)², φ(s) = 0
only at s = 0, and ∫₀^s φ → ∞ as |s| → ∞, by the same proof with
V = Σ_i ∫₀^{x_i} φ. (Global existence then follows because V is
non-increasing and radially unbounded, so solutions stay in compact
sublevel sets.)
