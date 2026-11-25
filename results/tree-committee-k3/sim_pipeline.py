"""Exploratory finite-size simulation of the three-stage construction (not part of the proof).

Needs numpy, scipy, cvxpy (with the Clarabel solver); not in requirements.txt.
    python sim_pipeline.py N alpha [seed]
Stage parameters are chosen from the limit formulas in check.py (float version below),
backed off slightly from the optimum.
"""
import sys
import numpy as np
import cvxpy as cp
from scipy.stats import norm
from scipy.optimize import brentq

Phi, phi = norm.cdf, norm.pdf


def I2(lo, hi, c):
    P = Phi(hi) - Phi(lo); M = phi(lo) - phi(hi); S = P - (hi * phi(hi) - lo * phi(lo))
    return c * c * P - 2 * c * M + S


def soft(x, rho, k):
    A, B = k - x, k - x * rho
    c = x * x * Phi(A) + I2(A, B, k) + (x * rho) ** 2 * (1 - Phi(B))
    return c, c - (Phi(B) - Phi(A)), Phi(A)


def hard(x, rho, k):
    B = k - x * rho
    c = k * k * Phi(B) + 2 * k * phi(B) + Phi(B) - B * phi(B) + (x * rho) ** 2 * (1 - Phi(B))
    return c, c - Phi(B)


def best_rho(LH, LR, k, margin=2e-3):
    best = None
    for rho in np.linspace(0.12, 0.6, 241):
        f = lambda x: LR * soft(x, rho, k)[0] + (LH * hard(x, rho, k)[0] if LH else 0) - 1
        x = brentq(f, 1e-4, 20)
        P = LR * soft(x, rho, k)[1] + (LH * hard(x, rho, k)[1] if LH else 0)
        e = soft(x, rho, k)[2]
        if P > margin and (best is None or e < best[0]):
            best = (e, rho, x)
    return best


def solve(X, H, R, rho, k):
    w = cp.Variable(X.shape[1]); h = X @ w
    obj = -rho * cp.sum(h)
    cons = [cp.norm(w) <= 1]
    if len(R):
        obj += (1 - rho) * cp.sum(cp.pos(k - h[R]))
    if len(H):
        cons.append(h[H] >= k)
    cp.Problem(cp.Minimize(obj), cons).solve(solver=cp.CLARABEL)
    return X @ w.value


def main():
    N, alpha = int(sys.argv[1]), float(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    k = 0.02; beta = 3 * alpha; m = int(beta * N)
    rng = np.random.default_rng(seed)
    X = [rng.standard_normal((m, N)) for _ in range(3)]
    e1p, rho1, _ = best_rho(0, beta, k)
    h1 = solve(X[0], [], np.arange(m), rho1, k); E1 = np.where(h1 <= 0)[0]
    print(f"stage 1: error fraction {len(E1)/m:.4f} (limit bound {e1p:.4f})")
    R = np.setdiff1d(np.arange(m), E1); a = len(E1) / N
    e2p, rho2, _ = best_rho(a, len(R) / N, k)
    h2 = solve(X[1], E1, R, rho2, k); E2 = R[h2[R] <= 0]
    print(f"stage 2: hard set satisfied {bool(np.all(h2[E1] > 0))}, error fraction {len(E2)/len(R):.4f} (limit bound {e2p:.4f})")
    T = np.concatenate([E1, E2]); print(f"stage 3 load |T|/N = {len(T)/N:.4f} (Cover needs < 2)")
    w3 = cp.Variable(N); p = cp.Problem(cp.Minimize(0), [X[2][T] @ w3 >= 1]); p.solve(solver=cp.CLARABEL)
    print("stage 3:", p.status)
    if w3.value is not None:
        Hm = np.stack([h1, h2, X[2] @ w3.value])
        print("all patterns stored:", bool(np.all(np.median(Hm, axis=0) > 0)))


if __name__ == "__main__":
    main()
