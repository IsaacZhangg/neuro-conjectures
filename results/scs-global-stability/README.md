# Global stability of the SCS random network

In the Sompolinsky–Crisanti–Sommers network dx/dt = −x + J tanh(x) with
J_ij ~ N(0, g²/N), the silent state x = 0 attracts every initial condition
with high probability for all g < 1/√2 ≈ 0.707, past the g < 1/2 that the
spectral-norm contraction argument gives. The conjectured range is g < 1.

Status: Partial
Kind: Theorem
Reproducibility: 🟡 written proof with runnable numerical spot checks; not yet independently reviewed

## Result

**Proved** (PROOF.md, Theorems 1 and 2): for every N and every s > 0, with
probability at least 1 − exp(−N s²/2), the origin is globally asymptotically
stable simultaneously for every g < 1/(√2 + s). For a fixed g < 1/√2 this
probability is at least 1 − exp(−N(1 − √2 g)²/(2g²)). Convergence is
exponential, at a rate that depends on g and the starting point.

The proof uses the Lyapunov function Σ log cosh x_i, which is the identity-
weight case of a classical criterion (Matsuoka 1992; Forti and Tesi 1995),
together with a non-asymptotic bound on the top eigenvalue of the symmetric
part of J. We found no paper that states this range for the SCS network.

**Still open:** 1/√2 ≤ g < 1. The identity-weight Lyapunov function cannot go
further, because the symmetric part's top eigenvalue tends to √2 g. In the
spot checks, the best diagonal weight found narrows the gap at finite N (from
g ≈ 0.706 to 0.738 at N = 1000), but the gain shrinks as N grows, so it may
vanish in the limit.

**Numerical spot checks** (check-output.json, guide only): the top eigenvalue
of (G + Gᵀ)/(2√N) averaged 1.389, 1.401, 1.406 and 1.409 at N = 250 to 2000,
below √2, with no sample above √2 + 0.1. The proved Lyapunov inequality held at
all 600 test states. Simulations at N = 500 from wide initial states converged
to 0 for every g tested from 0.6 to 0.95, consistent with the conjecture.

## Files

- [CLAIM.md](CLAIM.md): exact statement and decision rule.
- [PROOF.md](PROOF.md): the proof, the limits of the method, what remains open.
- [SOURCES.md](SOURCES.md): primary source, results used, later work, searches.
- [METHODS.md](METHODS.md): proof strategy and the numerical spot checks.
- [check.py](check.py): spot checks. `python check.py` (about 4 minutes) or
  `python check.py --quick` (a few seconds).
- [check-output.json](check-output.json): recorded full run.
- [REVIEW.md](REVIEW.md): review record.

## Scope and limits

Gaussian couplings with variance g²/N, tanh activation, no inputs or noise.
The theorem covers every initial condition, not just typical ones. It says
nothing about 1/√2 ≤ g < 1, and the simulations there are evidence, not proof.
A model result only: it says nothing about real brains.
