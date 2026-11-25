# Methods

## Model

dx/dt = −x + J tanh(x), J = (g/√N) G, G iid N(0,1), self-couplings included.

## Proof strategy

Lyapunov function V(x) = Σ log cosh x_i (Theorem 1, deterministic), plus a
non-asymptotic bound on the top eigenvalue of the symmetric part of J
(Sudakov–Fernique for the mean, Gaussian concentration for the tail;
Theorem 2). Full argument in PROOF.md.

## Numerical spot checks (check.py)

These test the inequalities the proof uses and map what happens beyond the
proved range. They do not set the status.

1. λ_max of W = (G + Gᵀ)/(2√N) for N = 250, 500, 1000, 2000, 200 samples each:
   sample mean against √2 and the empirical tail at s = 0.1 against
   exp(−N s²/2).
2. dV/dt against −(1 − λ)|tanh x|² at 600 states of mixed type (Gaussian at
   scales 10⁻³ to 10^1.5, sparse large coordinates, aligned with the top
   eigenvector of the symmetric part), N = 1000, g = 0.69.
3. RK4 simulations (dt = 0.1, T = 400) from Gaussian initial states of RMS 3
   and 0.3, N = 500, g ∈ {0.6, 0.69, 0.75, 0.85, 0.9, 0.95}, 6 initial states
   each. Reports the largest final RMS.
4. Best diagonal Lyapunov weight: minimise λ_max of sym(D^{1/2} J D^{−1/2})
   over positive diagonal D at g = 1 (L-BFGS on a soft-max of eigenvalues
   with increasing sharpness), N = 125, 250, 500, 1000, two matrices each.
   The value reported is the exact λ_max at the D found, so it is a valid
   upper bound on the optimum for that matrix.

Exit code 0 when checks 1 and 2 hold on every sample. Seeds fixed (default 0;
check 4 uses seeds 1000·N + s).

## Robustness

Theorem 2's bound does not depend on the diagonal of J (PROOF.md,
Assumptions). The proof does not use any property of tanh beyond
s·tanh s ≥ tanh² s for the GAS part.

## Deviations from plan

None. The claim and decision rule were committed first (CLAIM.md).
