"""Task 2: Lur'e-Postnikov certificate  V = x^T M x + 2 sum_i d_i log cosh x_i,  phi = tanh(x).
Vdot = 2 x^T M(-x + J phi) + 2 phi^T D(-x + J phi) <= Vdot + 2 S(x,phi),
S = sum of multiplier terms that are >= 0 for phi = tanh:
  sector   : tau_i phi_i (x_i - phi_i)                       (Q = diag(tau), tau >= 0)
  repeated : + lam_ij (phi_i -/+ phi_j)((x_i -/+ x_j) - (phi_i -/+ phi_j))  (tanh odd, slope in [0,1],
             identical in every coordinate)  <=>  Q symmetric, diagonally dominant, Q_ii >= sum_j |Q_ij|.
Quadratic form of Vdot + 2S in z = (x, phi):
  F = [[-2M, MJ - D + Q], [J^T M - D + Q, DJ + J^T D - 2Q]]  required  F <= -s blkdiag(0, I),  s > 0,
  M >= 0, d >= 0, normalisation tr M + sum d = N.   (V >= 0, Vdot <= -s|tanh x|^2 + boundedness => GAS.)
Also: omega = 0 necessary condition  Q(I - J) + (I - J)^T Q > 0  (x = J phi kills the M and D terms).

usage: python lurie.py N seed
"""
import sys, time
import numpy as np
import cvxpy as cp
from common import *


def lp_margin(J, variant, solver="CLARABEL"):
    N = J.shape[0]
    M = cp.Variable((N, N), symmetric=True)
    d = cp.Variable(N)
    s = cp.Variable()
    cons = [M >> 0, d >= 0, cp.trace(M) + cp.sum(d) == N]
    if variant == "noM":
        cons.append(M == 0)
    if variant in ("sector", "noM"):
        tau = cp.Variable(N)
        Q = cp.diag(tau)
        cons.append(tau >= 0)
    else:  # repeated (diagonally dominant Q), optionally no integral terms
        Q = cp.Variable((N, N), symmetric=True)
        U = cp.Variable((N, N), symmetric=True)
        off = np.ones((N, N)) - np.eye(N)
        cons += [U >= Q, U >= -Q, cp.diag(Q) >= cp.sum(cp.multiply(off, U), axis=1)]
        if variant == "repeated_d0":
            cons.append(d == 0)
    D = cp.diag(d)
    B = M @ J - D + Q
    C = D @ J + J.T @ D - 2 * Q
    F = cp.bmat([[-2 * M, B], [B.T, (C + C.T) / 2]])
    E = np.zeros((2 * N, 2 * N)); E[N:, N:] = np.eye(N)
    cons.append((F + F.T) / 2 << -s * E)
    prob = cp.Problem(cp.Maximize(s), cons)
    prob.solve(solver=solver)
    sol = dict(M=M.value, d=d.value, Q=Q.value if hasattr(Q, "value") else None)
    return prob.status, s.value, sol


def omega0_margin(J, dd=True):
    N = J.shape[0]
    Q = cp.Variable((N, N), symmetric=True)
    s = cp.Variable()
    cons = [cp.trace(Q) == N]
    if dd:
        U = cp.Variable((N, N), symmetric=True)
        off = np.ones((N, N)) - np.eye(N)
        cons += [U >= Q, U >= -Q, cp.diag(Q) >= cp.sum(cp.multiply(off, U), axis=1)]
    else:
        cons += [Q == cp.diag(cp.diag(Q)), cp.diag(Q) >= 0]
    A = np.eye(N) - J
    X = Q @ A + A.T @ Q
    cons.append((X + X.T) / 2 >> s * np.eye(N))
    prob = cp.Problem(cp.Maximize(s), cons)
    prob.solve(solver="CLARABEL")
    return prob.status, s.value


def bisect(fun, lo, hi, tol=1e-3):
    hist = []
    while hi - lo > tol:
        g = (lo + hi) / 2
        st, sv = fun(g)[:2]
        ok = sv is not None and sv > 1e-7
        hist.append((g, st, None if sv is None else float(sv)))
        if ok:
            lo = g
        else:
            hi = g
    return lo, hi, hist


def verify(J, sol, ntest=20000, rng=None):
    """Sample check of Vdot <= 0 along the true vector field (x random, multiple scales)."""
    rng = rng or np.random.default_rng(1)
    N = J.shape[0]
    M, d = sol["M"], sol["d"]
    worst = -np.inf
    for sc in (0.1, 1.0, 3.0, 10.0):
        X = sc * rng.standard_normal((ntest // 4, N))
        Ph = np.tanh(X)
        F = -X + Ph @ J.T
        Vd = 2 * np.sum((X @ M) * F, 1) + 2 * np.sum(Ph * d[None] * F, 1)
        worst = max(worst, np.max(Vd / np.sum(Ph ** 2, 1)))
    return float(worst)


if __name__ == "__main__":
    N, seed = int(sys.argv[1]), int(sys.argv[2])
    variants = sys.argv[3].split(",") if len(sys.argv) > 3 else ["sector", "repeated", "repeated_d0"]
    J0 = J0_of(N, seed)
    t0 = time.time()
    a = float(np.max(np.linalg.eigvals(J0).real))
    hi0 = min(2.0, 0.999 / a)
    out = dict(N=N, seed=seed, alpha_J0=a)
    lo, hi, _ = bisect(lambda g: omega0_margin(g * J0, dd=False), 0.4, hi0, tol=1e-4)
    out["g_LDS"] = lo
    lo, hi, _ = bisect(lambda g: omega0_margin(g * J0, dd=True), 0.4, hi0, tol=1e-3)
    out["g_omega0_DD"] = [lo, hi]
    print(out, flush=True)
    for v in variants:
        t = time.time()
        lo, hi, hist = bisect(lambda g: lp_margin(g * J0, v), 0.4, hi0, tol=1e-3)
        st, sv, sol = lp_margin(lo * J0, v)
        res = dict(g=[lo, hi], status_at_lo=st, s_at_lo=float(sv), statuses=sorted(set(h[1] for h in hist)),
                   minEigM=float(np.linalg.eigvalsh(sol["M"])[0]), trM=float(np.trace(sol["M"])),
                   min_d=float(np.min(sol["d"])), sum_d=float(np.sum(sol["d"])),
                   sampled_max_Vdot_over_phi2=verify(lo * J0, sol), time=time.time() - t)
        out[v] = res
        print(v, res, flush=True)
    out["time"] = time.time() - t0
    save_json(f"lurie_N{N}_s{seed}.json", out)
    print("RESULT", out)
