"""Task 3: counterexample search for sector-only arguments.
(a) constant diagonal C in [0,1]^N: a(J0) = max_C spectral abscissa(J0 C)  (homogeneous, so
    JC - I unstable for some C  <=>  g > 1/a(J0)).  Projected gradient ascent with restarts.
(b) switching C(t) (vertices, piecewise constant on a grid): lower bound on the worst-case
    Lyapunov exponent mu(g) = (1/T) log rho( prod_k expm((g J0 Lam_k - I) dt) ) via a
    Pontryagin-type forward/backward iteration (bang-bang update of Lam_k from the sign of
    (J^T p_k)_j x_{k,j}).
usage: python sector_search.py const N1,N2,.. nseeds   |   python sector_search.py switch N seed g1,g2,..
"""
import sys, time
import numpy as np
from scipy.linalg import expm
from common import *


def abscissa_and_grad(J0, c):
    A = J0 * c[None, :]
    w, VR = np.linalg.eig(A)
    k = np.argmax(w.real)
    lam, v = w[k], VR[:, k]
    # left eigenvector
    wl, VL = np.linalg.eig(A.T)
    kl = np.argmin(np.abs(wl - lam))
    u = VL[:, kl]  # u^T A = lam u^T
    grad = ((u @ J0) * v / (u @ v)).real
    return lam.real, grad, lam


def max_abscissa(J0, rng, n_restart=6, iters=400, lr=0.05):
    N = J0.shape[0]
    best, bestc, bestlam = -np.inf, None, None
    starts = [np.ones(N)] + [rng.random(N) for _ in range(n_restart // 2)] + \
             [(rng.random(N) < 0.5).astype(float) * 0.9 + 0.05 for _ in range(n_restart - n_restart // 2)]
    for c in starts:
        c = c.copy()
        step = lr
        f, gr, lam = abscissa_and_grad(J0, c)
        for it in range(iters):
            cn = np.clip(c + step * gr / (np.linalg.norm(gr) + 1e-12) * np.sqrt(N) * 0.1, 0, 1)
            fn, grn, lamn = abscissa_and_grad(J0, cn)
            if fn >= f:
                c, f, gr, lam = cn, fn, grn, lamn
                step = min(step * 1.2, 1.0)
            else:
                step *= 0.5
                if step < 1e-5:
                    break
        if f > best:
            best, bestc, bestlam = f, c.copy(), lam
    return best, bestc, bestlam


def switching_rate(J0, g, Lams, dt):
    N = J0.shape[0]
    Phis = [expm((g * J0 * L[None, :] - np.eye(N)) * dt) for L in Lams]
    Phi = np.eye(N)
    for P in Phis:
        Phi = P @ Phi
    w = np.linalg.eigvals(Phi)
    return np.log(np.max(np.abs(w))) / (len(Lams) * dt), Phis, Phi


def worst_switching(J0, g, rng, T=10.0, dt=0.1, iters=40, init=None, restarts=3):
    N = J0.shape[0]
    K = int(round(T / dt))
    best = -np.inf
    bestL = None
    inits = []
    if init is not None:
        inits.append([init.copy() for _ in range(K)])
    inits.append([np.ones(N) for _ in range(K)])
    for _ in range(restarts):
        inits.append([(rng.random(N) < 0.5).astype(float) for _ in range(K)])
    for Lams in inits:
        for it in range(iters):
            mu, Phis, Phi = switching_rate(J0, g, Lams, dt)
            if mu > best:
                best, bestL = mu, [L.copy() for L in Lams]
            w, VR = np.linalg.eig(Phi)
            k = np.argmax(np.abs(w))
            v = VR[:, k]
            wl, VL = np.linalg.eig(Phi.T)
            u = VL[:, np.argmin(np.abs(wl - w[k]))]
            # forward states at interval midpoints (approx: start of interval), backward costates
            xs = [v]
            for P in Phis[:-1]:
                xs.append(P @ xs[-1])
            ps = [None] * K
            p = u
            for kk in range(K - 1, -1, -1):
                ps[kk] = p  # costate at end of interval kk
                p = Phis[kk].T @ p
            new = []
            for kk in range(K):
                xm = Phis[kk] @ xs[kk]  # state at end of interval
                xmid = 0.5 * (xs[kk] + xm)
                score = ((g * J0.T @ ps[kk]) * xmid / (u @ v)).real
                new.append((score > 0).astype(float))
            # damped update: change a random half of intervals
            flip = rng.random(K) < 0.5
            Lams = [new[kk] if flip[kk] else Lams[kk] for kk in range(K)]
    return best, bestL


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "const":
        Ns = [int(x) for x in sys.argv[2].split(",")]
        ns = int(sys.argv[3])
        out = {}
        for N in Ns:
            rows = []
            for seed in range(ns):
                J0 = J0_of(N, seed)
                rng = np.random.default_rng(seed)
                t = time.time()
                a1 = spectral_abscissa(J0)
                a, c, lam = max_abscissa(J0, rng)
                frac01 = float(np.mean((c < 1e-6) | (c > 1 - 1e-6)))
                rows.append(dict(seed=seed, alpha_C_eq_I=float(a1), a_max=float(a), g_C=float(1 / a),
                                 lam_imag=float(abs(lam.imag)), frac_c_at_bounds=frac01,
                                 mean_c=float(c.mean()), t=time.time() - t))
                print(N, rows[-1], flush=True)
                np.save(os.path.join(LOGDIR, f"Cstar_N{N}_s{seed}.npy"), c)
            out[N] = rows
        save_json("sector_const_" + "_".join(map(str, Ns)) + ".json", out)
    else:
        N, seed = int(sys.argv[2]), int(sys.argv[3])
        gs = [float(x) for x in sys.argv[4].split(",")]
        J0 = J0_of(N, seed)
        rng = np.random.default_rng(seed)
        cpath = os.path.join(LOGDIR, f"Cstar_N{N}_s{seed}.npy")
        c = np.load(cpath) if os.path.exists(cpath) else None
        out = dict(N=N, seed=seed, rows=[])
        for g in gs:
            t = time.time()
            const = (g * spectral_abscissa(J0 * c[None, :]) - 1) if c is not None else None
            mu, L = worst_switching(J0, g, rng, init=c)
            nsw = int(sum(np.sum(L[k] != L[k - 1]) for k in range(1, len(L))))
            row = dict(g=g, mu_switch=float(mu), mu_constCstar=const, mu_constI=float(g * spectral_abscissa(J0) - 1),
                       n_bitflips=nsw, t=time.time() - t)
            print(row, flush=True)
            out["rows"].append(row)
        save_json(f"sector_switch_N{N}_s{seed}.json", out)
