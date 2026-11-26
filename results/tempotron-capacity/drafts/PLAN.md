# Work in progress (paused 2026-11-26 at the coordinator's request)

Nothing here is a claim yet. CLAIM.md and METHODS.md are not frozen.

## Model
x^mu(t) iid N(0, I_N), t = 1..T, labels iid fair ±1, free (w, theta);
output + iff max_t w·x(t) > theta. alpha_c(T) as in NOVELTY.md.

## Upper bound A (folklore, record only)
Cover counting over the PT vectors: alpha_c(T) <= alpha* with
T·H2(1/(alpha T)) = 1; equals 2 at T = 1, about log2 T + log2 log2 T + log2 e.

## Upper bound B (planned, new)
epsilon-net over w (eps ~ c/sqrt(ln T)) and theta (spacing ~ c2/sqrt(ln T));
uniform operator-norm bound over m-subsets limits patterns an adaptive change
can flip to rho N; other patterns flip only if max_t w0·x(t) is within eta of
theta0 (probability ~ eta·density of max of T Gaussians); union bound gives
alpha_c(T) <= (1/2) log2 ln T + C, matching the leading term of Rubin,
Monasson and Sompolinsky (2010). Constants still to be pinned.

## Lower bound (planned, new)
Random direction u independent of data; t_mu = argmax_t u·x^mu(t) for each
positive pattern; w = a u + v with v ⊥ u. Conditional on the u-projections the
orthogonal parts are iid Gaussian, so storage is a random linear feasibility
problem in v. CGMT (Thrampoulidis, Oymak and Hassibi 2015), with the
lambda-ball radius L -> infinity trick to get exact feasibility, gives:
storable w.h.p. if some A, Theta satisfy
  alpha·(1/2)·[E(Theta − A s − z)_+^2 + T·E(sigma z − Theta)_+^2] < 1,
s = max of T iid N(0,1), z ~ N(0,1) independent, sigma^2 = 1 + A^2.
Gives exactly 2 at T = 1. drafts/lb_explore.py evaluates it for larger T
(run was stopped by the pause before T >= 2 finished; it is slow, the quad
call should be vectorised).

## Next steps when resumed
1. Finish lb_explore numerics (fast version); see whether the bound exceeds 2
   and how it grows with T.
2. Freeze CLAIM.md and METHODS.md.
3. PROOF.md for both bounds, certified check.py (python-flint), --quick mode.
4. One independent agent review, catalogue checks, README/CONTENTS/HISTORY,
   merge to main.
