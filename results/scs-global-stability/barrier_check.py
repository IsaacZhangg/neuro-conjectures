"""Numerical guide for BARRIER.md (Theorem 3). Not part of the proof.

For each (N, g) it draws one Gaussian matrix G (seed 2026), sets J = g G / sqrt(N),
A = J - I, and tests two positive semidefinite certificates X. If every diagonal
entry of A X is >= 0, then J - I is not Lyapunov diagonally stable (Lemma 4).
Run: python barrier_check.py   (writes barrier-check-output.json)
"""
import json
import numpy as np

SIZES = [500, 1000, 2000, 4000]
GAINS = [0.8, 0.9, 1.0]


def summarize(d):
    return dict(min_over_mean=float(d.min() / d.mean()), frac_negative=float((d < 0).mean()),
                certifies_not_LDS=bool(d.min() >= 0))


def main():
    rng = np.random.default_rng(2026)
    rows = []
    for N in SIZES:
        for g in GAINS:
            G = rng.standard_normal((N, N))
            J = g * G / np.sqrt(N)
            A = J - np.eye(N)
            Js = (J + J.T) / 2
            K = (J - J.T) / 2
            lam, V = np.linalg.eigh(Js)
            a = (1 + np.sqrt(2) * g) / 2
            sel = lam >= a
            m = int(sel.sum())
            row = dict(N=N, g=g, a=float(a), m=m, lambda_max=float(lam[-1]))
            if m:
                P = V[:, sel] @ V[:, sel].T
                dP = np.einsum("ij,ji->i", A, P)
                T = np.einsum("ij,ji->i", K, P)
                Pii = np.diag(P)
                row["projection"] = summarize(dP)
                row["projection"]["min_Pii_times_N_over_m"] = float(Pii.min() * N / m)
                row["projection"]["max_T_over_margin"] = float(np.max(np.abs(T) / ((a - 1) * Pii)))
            w = np.clip(lam - 1, 0, None)
            X = (V * w) @ V.T
            row["weighted"] = summarize(np.einsum("ij,ji->i", A, X))
            row["weighted"]["rank"] = int((w > 0).sum())
            rows.append(row)
            print(json.dumps(row))
    with open("barrier-check-output.json", "w") as fh:
        json.dump(dict(seed=2026, rows=rows), fh, indent=1)


if __name__ == "__main__":
    main()
