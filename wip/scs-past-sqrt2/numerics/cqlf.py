"""Task 1: common quadratic Lyapunov function (Krasovskii / Euclidean-P contraction) for
A(Lambda) = -I + g J0 Lambda, Lambda diagonal in [0,1]^N (vertex check on Lambda = diag(1_S)).
Cutting-plane SDP inside a bisection on g.  Also: best diagonal P (LDS of J - I) threshold,
diagonal-P vertex CQLF threshold, and P = I threshold.

usage: python cqlf.py N seed [mode]   mode in {full, diag}  (default full)
"""
import sys, time
import numpy as np
import cvxpy as cp
from common import *


def solve_ws(J0, g, W, diag=False, solver="CLARABEL"):
    N = J0.shape[0]
    if diag:
        d = cp.Variable(N)
        P = cp.diag(d)
    else:
        P = cp.Variable((N, N), symmetric=True)
    s = cp.Variable()
    cons = [cp.trace(P) == N]
    I = np.eye(N)
    for m in W:
        B = g * J0 * m[None, :]
        X = -2 * P + P @ B + B.T @ P
        cons.append((X + X.T) / 2 << -s * I)
    prob = cp.Problem(cp.Maximize(s), cons)
    try:
        prob.solve(solver=solver)
    except Exception as e:
        prob.solve(solver="SCS", eps=1e-7, max_iters=20000)
    Pv = P.value if not diag else np.diag(d.value)
    return prob.status, (s.value if s.value is not None else -np.inf), (None if Pv is None else (Pv + Pv.T) / 2)


def feasible_at(J0, g, W, rng, diag=False, max_rounds=200, exhaustive=False, log=None):
    """Cutting-plane loop. Returns (feasible:bool, P, info)."""
    N = J0.shape[0]
    info = dict(rounds=0, heur_missed=0)
    for r in range(max_rounds):
        info["rounds"] = r + 1
        st, sv, P = solve_ws(J0, g, W, diag=diag)
        info["status"] = st
        if P is None or sv is None or sv < 1e-7:
            info["s"] = float(sv) if sv is not None else None
            return False, P, info
        info["s"] = float(sv)
        cands = heuristic_max(P, J0, g, extra_starts=W[-5:], n_random=20 if N > 40 else 30, rng=rng)
        viol = [m for f, m in cands if f > -0.25 * sv][:10]
        worst = cands[0][0]
        if worst >= 0:
            W.extend(viol if viol else [cands[0][1]])
            continue
        if exhaustive:
            ex, exm = exhaustive_max(P, J0, g)
            info["exh_max"] = float(ex)
            info["heur_max"] = float(worst)
            if ex >= 0:
                info["heur_missed"] += 1
                W.append(exm)
                continue
            if ex > worst + 1e-9:
                info["heur_gap"] = float(ex - worst)
        info["cert_margin"] = float(-worst)
        return True, P, info
    return None, P, info


def bisect(J0, lo, hi, rng, diag=False, tol=2e-3, exhaustive=False, log=print, tcap=1e9):
    N = J0.shape[0]
    W = [np.zeros(N), np.ones(N)] + [np.eye(N)[i] for i in range(N)]
    # top-eigvec sign patterns of sym(J0)
    w, V = np.linalg.eigh((J0 + J0.T) / 2)
    for k in range(1, 4):
        W.append((V[:, -k] * (J0.T @ V[:, -k]) > 0).astype(float))
    ok, P, info = feasible_at(J0, lo, W, rng, diag=diag, exhaustive=exhaustive)
    log(f"  lo={lo:.4f} feasible={ok} {info}")
    if not ok:
        return None, lo, W, P
    Pbest = P
    tstart = time.time()
    while hi - lo > tol and time.time() - tstart < tcap:
        g = (lo + hi) / 2
        t = time.time()
        ok, P, info = feasible_at(J0, g, W, rng, diag=diag, exhaustive=exhaustive)
        log(f"  g={g:.4f} feasible={ok} |W|={len(W)} t={time.time()-t:.1f}s {info}")
        if ok:
            lo, Pbest = g, P
        else:
            hi = g
    return lo, hi, W, Pbest


def lds_threshold(J0, tol=1e-4):
    """Largest g with diagonal D>0: D(gJ0 - I) + (.)^T D < 0  (single LMI, Lambda = I)."""
    N = J0.shape[0]
    lo, hi = 0.3, 1.0 / max(np.max(np.linalg.eigvals(J0).real), 1e-9)
    hi = min(hi, 1.2)
    while hi - lo > tol:
        g = (lo + hi) / 2
        st, sv, P = solve_ws(J0, g, [np.ones(N), np.zeros(N)], diag=True)
        if sv is not None and sv > 1e-8:
            lo = g
        else:
            hi = g
    return lo, hi


def PI_threshold(J0, rng):
    """P = I: g < 1 / max_S lambda_max(sym(J0 Lambda_S))."""
    N = J0.shape[0]
    P = np.eye(N)
    if N <= 16:
        v, m = exhaustive_max(P, J0, 1.0)
    else:
        v, m = heuristic_max(P, J0, 1.0, n_random=50, rng=rng)[0]
    # L = -2I + 2 sym(J0 Lam) -> lmax = -2 + 2 mu
    mu = (v + 2) / 2
    return 1.0 / mu, mu


if __name__ == "__main__":
    N, seed = int(sys.argv[1]), int(sys.argv[2])
    mode = sys.argv[3] if len(sys.argv) > 3 else "full"
    J0 = J0_of(N, seed)
    rng = np.random.default_rng(seed)
    out = dict(N=N, seed=seed, mode=mode)
    t0 = time.time()
    lines = []
    def log(s):
        print(s, flush=True); lines.append(s)
    out["alpha_J0"] = float(np.max(np.linalg.eigvals(J0).real))
    out["lmax_symJ0"] = float(np.linalg.eigvalsh((J0 + J0.T) / 2)[-1])
    out["g_LDS"] = lds_threshold(J0)
    out["g_PI"] = PI_threshold(J0, rng)[0]
    log(f"N={N} seed={seed} alpha(J0)={out['alpha_J0']:.4f} lmax sym={out['lmax_symJ0']:.4f} "
        f"g_LDS={out['g_LDS'][0]:.4f} g_P=I={out['g_PI']:.4f}")
    diag = (mode == "diag")
    hi0 = min(1.0, 0.999 / out["alpha_J0"])
    gl = out["g_LDS"][0]
    tcap = float(sys.argv[4]) if len(sys.argv) > 4 else 600
    lo, hi, W, P = bisect(J0, gl - 0.02, min(gl + 0.02, hi0), rng, diag=diag, tol=2.5e-3,
                          exhaustive=(N <= 16), log=log, tcap=tcap)
    if lo is None:
        log("  lower bracket infeasible -> restart bracket from 0.45")
        lo, hi, W, P = bisect(J0, 0.45, gl - 0.02, rng, diag=diag, tol=2.5e-3,
                              exhaustive=(N <= 16), log=log, tcap=tcap)
    out["g_cqlf"] = [lo, hi]
    out["W_size"] = len(W)
    if P is not None and lo is not None:
        ev = np.linalg.eigvalsh(P)
        out["P_cond"] = float(ev[-1] / ev[0])
        offd = P - np.diag(np.diag(P))
        out["P_offdiag_frac"] = float(np.linalg.norm(offd) / np.linalg.norm(P))
        # final validation of certificate at g=lo
        if N <= 16:
            ex, _ = exhaustive_max(P, J0, lo)
            out["final_exhaustive_max"] = float(ex)
        else:
            out["final_heur_max"] = float(heuristic_max(P, J0, lo, n_random=200, rng=np.random.default_rng(99))[0][0])
        np.save(os.path.join(LOGDIR, f"P_{mode}_N{N}_s{seed}.npy"), P)
    out["time"] = time.time() - t0
    log(f"RESULT {out}")
    save_json(f"cqlf_{mode}_N{N}_s{seed}.json", out)
