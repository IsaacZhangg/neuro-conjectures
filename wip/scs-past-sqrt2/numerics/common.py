"""Shared helpers for the SCS certificate numerics.
Model: dx/dt = -x + J tanh(x), J = g * J0, J0 = G/sqrt(N), G iid N(0,1).
"""
import numpy as np
import itertools, time, json, os

LOGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
os.makedirs(LOGDIR, exist_ok=True)


def J0_of(N, seed):
    rng = np.random.default_rng(1000 * N + seed)
    return rng.standard_normal((N, N)) / np.sqrt(N)


def L_S(P, J0, g, mask):
    """P A + A^T P with A = -I + g J0 diag(mask)."""
    K = g * (P @ J0) * mask[None, :]
    return -2 * P + K + K.T


def lmax(M):
    return np.linalg.eigvalsh(M)[-1]


def all_masks(N):
    # (2^N, N) array of 0/1
    idx = np.arange(2 ** N, dtype=np.int64)[:, None]
    return ((idx >> np.arange(N)) & 1).astype(float)


def exhaustive_max(P, J0, g, batch=8192):
    """max over all 2^N vertex masks of lambda_max(L_S(P)). Returns (value, mask)."""
    N = P.shape[0]
    PJ = g * (P @ J0)
    best, bestm = -np.inf, None
    for start in range(0, 2 ** N, batch):
        idx = np.arange(start, min(2 ** N, start + batch), dtype=np.int64)[:, None]
        m = ((idx >> np.arange(N)) & 1).astype(float)
        K = PJ[None, :, :] * m[:, None, :]
        L = -2 * P[None] + K + np.transpose(K, (0, 2, 1))
        ev = np.linalg.eigvalsh(L)[:, -1]
        k = np.argmax(ev)
        if ev[k] > best:
            best, bestm = ev[k], m[k].copy()
    return best, bestm


def local_search(P, J0, g, mask0, max_iter=200):
    """Ascent on f(S) = lambda_max(L_S(P)): alternate (eigvec -> best S for it),
    then single-flip hill climbing with exact eigenvalues."""
    N = P.shape[0]
    PJ = g * (P @ J0)
    m = mask0.copy()
    f = None
    for _ in range(max_iter):
        # alternating step
        for _a in range(50):
            w, V = np.linalg.eigh(L_S(P, J0, g, m))
            v = V[:, -1]
            f_new = w[-1]
            score = (PJ.T @ v) * v  # coefficient of mask_j in v^T L v / 2
            m_new = (score > 0).astype(float)
            if np.array_equal(m_new, m):
                break
            m = m_new
        f = lmax(L_S(P, J0, g, m))
        # single flips (batched exact)
        flips = np.repeat(m[None], N, 0)
        flips[np.arange(N), np.arange(N)] = 1 - flips[np.arange(N), np.arange(N)]
        K = PJ[None] * flips[:, None, :]
        L = -2 * P[None] + K + np.transpose(K, (0, 2, 1))
        ev = np.linalg.eigvalsh(L)[:, -1]
        k = np.argmax(ev)
        if ev[k] > f + 1e-12:
            m = flips[k].copy()
        else:
            break
    return f, m


def heuristic_max(P, J0, g, extra_starts=(), n_random=30, rng=None):
    """Heuristic search for the most violated vertex. Returns list of (value, mask) sorted desc."""
    N = P.shape[0]
    rng = rng or np.random.default_rng(0)
    starts = [np.ones(N), np.zeros(N)]
    # sign patterns of top eigvecs of L_I(P)
    w, V = np.linalg.eigh(L_S(P, J0, g, np.ones(N)))
    PJ = g * (P @ J0)
    for k in range(1, 4):
        v = V[:, -k]
        starts.append(((PJ.T @ v) * v > 0).astype(float))
    starts += [np.asarray(s, float) for s in extra_starts]
    for _ in range(n_random):
        starts.append((rng.random(N) < rng.uniform(0.2, 0.8)).astype(float))
    res = {}
    for s in starts:
        f, m = local_search(P, J0, g, s)
        res[tuple(m.astype(int))] = f
    out = sorted([(f, np.array(k, float)) for k, f in res.items()], key=lambda t: -t[0])
    return out


def spectral_abscissa(A):
    return np.max(np.linalg.eigvals(A).real)


def save_json(name, obj):
    with open(os.path.join(LOGDIR, name), "w") as f:
        json.dump(obj, f, indent=1, default=float)
