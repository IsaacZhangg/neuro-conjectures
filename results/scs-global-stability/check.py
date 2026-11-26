"""Numerical spot checks for the SCS global-stability result.

The result itself is the written proof in PROOF.md. This script only checks,
on random samples, the inequalities the proof uses and maps what happens
past the proved range. It never sets the status on its own.

Checks
  1. lambda_max of W = (G + G^T) / (2 sqrt N): sample mean <= sqrt 2
     (Sudakov-Fernique bound) and empirical tail vs exp(-N s^2 / 2).
  2. Lyapunov derivative dV/dt <= -(1 - lambda) |tanh x|^2 at random states
     (V = sum log cosh x_i), for a g just inside the proved range.
  3. Simulations from wide initial conditions for g in and past the proved
     range (guide only).
  4. Best diagonal Lyapunov weight D found by optimisation: how far the
     diagonal-weight method could reach at finite N (guide only).

Exit 0 when every proved inequality holds on every sample, 1 otherwise.
"""
import argparse
import json
import platform
import sys
import time

import numpy as np
import scipy
from scipy.optimize import minimize


def lam_sym(J):
    return np.linalg.eigvalsh((J + J.T) / 2)[-1]


def check_edge(rng, sizes, reps):
    out = []
    ok = True
    for N in sizes:
        vals = np.array([lam_sym(rng.standard_normal((N, N)) / np.sqrt(N)) for _ in range(reps)])
        # tail at s = 0.1 versus the proved bound exp(-N s^2 / 2)
        s = 0.1
        emp = float(np.mean(vals >= np.sqrt(2) + s))
        bound = float(np.exp(-N * s * s / 2))
        row = dict(N=N, reps=reps, mean=float(vals.mean()), sd=float(vals.std()),
                   max=float(vals.max()), tail_s0p1_empirical=emp, tail_s0p1_bound=bound)
        # the mean bound is E <= sqrt2; allow 3 standard errors of sampling noise
        if vals.mean() > np.sqrt(2) + 3 * vals.std() / np.sqrt(reps):
            ok = False
        if emp > bound + 3 * np.sqrt(max(bound, 1e-12) / reps) + 1.0 / reps:
            ok = False
        out.append(row)
    return ok, out


def check_lyapunov(rng, N, g, states):
    J = g * rng.standard_normal((N, N)) / np.sqrt(N)
    lam = lam_sym(J)
    worst = -np.inf
    for k in range(states):
        scale = 10 ** rng.uniform(-3, 1.5)
        x = scale * rng.standard_normal(N)
        if k % 3 == 1:  # sparse large coordinates
            x = x * (rng.random(N) < 0.05)
        if k % 3 == 2:  # aligned with the top eigenvector of the symmetric part
            w, V = np.linalg.eigh((J + J.T) / 2)
            x = np.arctanh(np.clip(scale / (1 + scale) * V[:, -1] / np.abs(V[:, -1]).max(), -0.999, 0.999))
        phi = np.tanh(x)
        dV = phi @ (-x + J @ phi)
        bound = -(1 - lam) * (phi @ phi)
        if phi @ phi > 0:
            worst = max(worst, (dV - bound) / (phi @ phi))
    return bool(worst <= 1e-10), dict(N=N, g=g, lambda_max_sym=float(lam), states=states,
                                       worst_excess_over_bound=float(worst))


def simulate(rng, N, g, inits, T, dt=0.1):
    J = g * rng.standard_normal((N, N)) / np.sqrt(N)
    f = lambda x: -x + J @ np.tanh(x)
    finals = []
    for k in range(inits):
        x = (3.0 if k % 2 == 0 else 0.3) * rng.standard_normal(N)
        for _ in range(int(T / dt)):
            k1 = f(x); k2 = f(x + dt / 2 * k1); k3 = f(x + dt / 2 * k2); k4 = f(x + dt * k3)
            x = x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        finals.append(float(np.sqrt(np.mean(x ** 2))))
    return dict(N=N, g=g, inits=inits, T=T, lambda_max_sym=float(lam_sym(J)),
                spectral_radius=float(np.max(np.abs(np.linalg.eigvals(J)))),
                max_final_rms=max(finals))


def best_diagonal(J, betas=(20, 60, 200, 600, 2000)):
    """Minimise lambda_max of sym(D^1/2 J D^-1/2) over positive diagonal D.

    Returns a valid upper bound on the minimum (the value at the D found)."""
    N = len(J)
    h = np.zeros(N)
    for beta in betas:
        def fun(h):
            h = h - h.mean()
            e = np.exp(h / 2)
            C = e[:, None] * J / e[None, :]
            w, V = np.linalg.eigh((C + C.T) / 2)
            m = w[-1]
            p = np.exp(beta * (w - m)); s = p.sum(); p /= s
            grad = 0.5 * (V * (C @ V - C.T @ V)) @ p
            return m + np.log(s) / beta, grad - grad.mean()
        h = minimize(fun, h, jac=True, method="L-BFGS-B", options={"maxiter": 6000}).x
    h = h - h.mean()
    e = np.exp(h / 2)
    C = e[:, None] * J / e[None, :]
    return float(np.linalg.eigvalsh((C + C.T) / 2)[-1]), float(h.std())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    q = args.quick
    ok1, edge = check_edge(rng, [100, 200] if q else [250, 500, 1000, 2000], 20 if q else 200)
    ok2, lyap = check_lyapunov(rng, 100 if q else 1000, 0.69, 60 if q else 600)
    sims = [simulate(rng, 100 if q else 500, g, 2 if q else 6, 60 if q else 400)
            for g in ([0.6, 0.9] if q else [0.6, 0.69, 0.75, 0.85, 0.9, 0.95])]
    diag = []
    for N in ([40] if q else [125, 250, 500, 1000]):
        for s in range(1 if q else 2):
            J = np.random.default_rng(1000 * N + s).standard_normal((N, N)) / np.sqrt(N)
            lD, hstd = best_diagonal(J, betas=(20, 200) if q else (20, 60, 200, 600, 2000))
            ls = lam_sym(J)
            diag.append(dict(N=N, seed=s, lambda_sym=float(ls), lambda_bestD=lD,
                             gap=float(ls - lD), g_sym=float(1 / ls), g_bestD=float(1 / lD),
                             log_d_sd=hstd))
    passed = ok1 and ok2
    out = dict(result="scs-global-stability", quick=q, seed=args.seed, passed=passed,
               proved_range="g < 1/sqrt(2) = %.6f, probability >= 1 - exp(-N(1 - sqrt2 g)^2/(2 g^2))" % (1 / np.sqrt(2)),
               edge_check=edge, lyapunov_check=lyap, simulations=sims, best_diagonal_weight=diag,
               versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
               runtime_s=round(time.time() - t0, 1))
    if not q:
        with open("check-output.json", "w") as fh:
            json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
