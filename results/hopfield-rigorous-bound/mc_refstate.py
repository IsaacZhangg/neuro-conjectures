"""Monte Carlo sanity check for PROOF-refstate.md (exploratory, not a proof).

For model G with xi = 1 (gauge), builds s = G1, s_hat, z = (GP)^T s_hat,
w = z + nu*1 with nu ~ N(0, 1/N), T = {i : w_i < -theta}, sigma1 = 1 - 2 e_T,
and checks
  (a) the exact identity  ||G sigma1||^2 = (||s||(1-2p) + 2N zeta)^2 + 4||R e_T||^2
      (floating point, relative error should be ~1e-12),
  (b) the heuristic limits p -> Q(theta), zeta -> phi(theta), ||s||^2/N^2 -> m,
      ||R e_T||^2/N^2 -> m p(1-p), and
      gain := [(1.sigma1)^2 + ||G sigma1||^2]/N^2 - (1 + ||s||^2/N^2) -> gamma(m, theta).
Usage: python mc_refstate.py [--N 2000] [--trials 5] [--m 0.13] [--theta 2.53]
"""
import argparse
import numpy as np
from scipy.stats import norm


def one_trial(rng, N, M, theta):
    G = rng.standard_normal((M, N))
    one = np.ones(N)
    s = G @ one
    s_norm = np.linalg.norm(s)
    s_hat = s / s_norm
    GP = G - np.outer(s, one) / N          # G P with P = I - 11^T/N
    z = GP.T @ s_hat
    R = GP - np.outer(s_hat, z)            # (I - s_hat s_hat^T) G P
    nu = rng.standard_normal() / np.sqrt(N)
    w = z + nu
    T = w < -theta
    eT = T.astype(float)
    p = eT.sum() / N
    zeta = -z[T].sum() / N
    sigma1 = one - 2 * eT
    lhs = np.linalg.norm(G @ sigma1) ** 2
    ReT = R @ eT
    rhs = (s_norm * (1 - 2 * p) + 2 * N * zeta) ** 2 + 4 * ReT @ ReT
    gain = ((N * (1 - 2 * p)) ** 2 + lhs) / N**2 - (1 + s_norm**2 / N**2)
    return {
        "rel_err_identity": abs(lhs - rhs) / lhs,
        "ReT_perp_shat": abs(ReT @ s_hat) / (np.linalg.norm(ReT) + 1e-300),
        "p": p, "zeta": zeta, "s2": s_norm**2 / N**2,
        "ReT2_over_mp1p": (ReT @ ReT) / N**2 / ((M / N) * p * (1 - p)),
        "gain": gain,
        "d_sigma1": int(T.sum()),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=2000)
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--m", type=float, default=0.13)
    ap.add_argument("--theta", type=float, default=2.53)
    ap.add_argument("--seed", type=int, default=20261125)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    N, M, th = a.N, int(round(a.m * a.N)), a.theta
    m = M / N
    Qt, ft = norm.sf(th), norm.pdf(th)
    gamma = 4 * (np.sqrt(m) * ft * (1 - 2 * Qt) + ft**2 - Qt * (1 - Qt))
    print(f"N={N} M={M} m={m:.4f} theta={th}")
    print(f"limits: p->{Qt:.6f} zeta->{ft:.6f} ||s||^2/N^2->{m:.6f} "
          f"||Re_T||^2/(N^2 m p(1-p))->1 gain->gamma={gamma:.6f}")
    rows = [one_trial(rng, N, M, th) for _ in range(a.trials)]
    for r in rows:
        print("  " + "  ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}"
                               for k, v in r.items()))
    for k in ("p", "zeta", "s2", "ReT2_over_mp1p", "gain"):
        v = np.array([r[k] for r in rows])
        print(f"mean {k:15s} = {v.mean():.6f}  (sd {v.std(ddof=1):.2g})")
    worst = max(r["rel_err_identity"] for r in rows)
    print(f"max relative error of exact identity: {worst:.2e}")
    return 0 if worst < 1e-9 else 1


if __name__ == "__main__":
    raise SystemExit(main())
