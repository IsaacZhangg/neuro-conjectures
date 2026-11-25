"""Tree measurements and the scaling-law fit used by check.py.

Everything here works on a tree given as vertex coordinates (micrometres) and
undirected edges. No data access happens in this module.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import ConvexHull, QhullError


@dataclass
class TreeStats:
    L: float            # total cable length after pruning, um
    n_branch: int       # nodes of degree >= 3 after pruning
    n_tips: int         # nodes of degree 1 after pruning
    V_hull: float       # convex hull volume of the pruned tree's vertices, um^3
    n_vertices: int
    n_components: int


def _degrees(n_vertices: int, edges: np.ndarray) -> np.ndarray:
    deg = np.zeros(n_vertices, dtype=np.int64)
    np.add.at(deg, edges[:, 0], 1)
    np.add.at(deg, edges[:, 1], 1)
    return deg


def prune_terminal_segments(xyz: np.ndarray, edges: np.ndarray, tau: float):
    """Remove, in a single pass, every terminal segment shorter than tau (um).

    A terminal segment is the path from a tip (degree 1) to the first vertex
    of degree >= 3. Isolated paths with no branch point are never removed.
    Returns the kept vertex coordinates and re-indexed edges.
    """
    n = len(xyz)
    if tau <= 0 or len(edges) == 0:
        return xyz, edges
    deg = _degrees(n, edges)
    seglen = np.linalg.norm(xyz[edges[:, 0]] - xyz[edges[:, 1]], axis=1)
    # adjacency lists in CSR form
    order = np.argsort(np.concatenate([edges[:, 0], edges[:, 1]]), kind="stable")
    nbr = np.concatenate([edges[:, 1], edges[:, 0]])[order]
    elen = np.concatenate([seglen, seglen])[order]
    start = np.zeros(n + 1, dtype=np.int64)
    np.add.at(start, np.concatenate([edges[:, 0], edges[:, 1]]) + 1, 1)
    start = np.cumsum(start)
    remove = np.zeros(n, dtype=bool)
    for tip in np.flatnonzero(deg == 1):
        path = [tip]
        prev, cur, length = -1, tip, 0.0
        while True:
            nb = nbr[start[cur]:start[cur + 1]]
            ln = elen[start[cur]:start[cur + 1]]
            nxt = [(v, l) for v, l in zip(nb, ln) if v != prev]
            if not nxt:          # reached another tip: unbranched path, keep
                path = None
                break
            v, l = nxt[0]
            length += l
            if length >= tau:
                path = None
                break
            if deg[v] >= 3:      # reached branch point: segment is short
                break
            path.append(v)
            prev, cur = cur, v
        if path is not None:
            remove[path] = True
    keep = ~remove
    new_index = -np.ones(n, dtype=np.int64)
    new_index[keep] = np.arange(keep.sum())
    e_keep = keep[edges[:, 0]] & keep[edges[:, 1]]
    return xyz[keep], new_index[edges[e_keep]]


def tree_stats(xyz: np.ndarray, edges: np.ndarray, tau: float) -> TreeStats:
    xyz, edges = prune_terminal_segments(np.asarray(xyz, float), np.asarray(edges, np.int64), tau)
    n = len(xyz)
    deg = _degrees(n, edges)
    L = float(np.linalg.norm(xyz[edges[:, 0]] - xyz[edges[:, 1]], axis=1).sum()) if len(edges) else 0.0
    try:
        V = float(ConvexHull(xyz).volume) if n >= 4 else 0.0
    except QhullError:
        V = 0.0
    n_comp = n - len(edges)  # forest: components = vertices - edges
    return TreeStats(L, int((deg >= 3).sum()), int((deg == 1).sum()), V, n, int(n_comp))


def fit_scaling(L, V, n):
    """OLS fit of log L = c + a log V + b log n. Returns (a, b, c)."""
    X = np.column_stack([np.ones(len(L)), np.log(V), np.log(n)])
    coef, *_ = np.linalg.lstsq(X, np.log(L), rcond=None)
    return float(coef[1]), float(coef[2]), float(coef[0])


def bootstrap_fit(L, V, n, strata, n_boot: int, rng: np.random.Generator):
    """Percentile 95% intervals for (a, b), resampling neurons within strata."""
    L, V, n, strata = map(np.asarray, (L, V, n, strata))
    groups = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    out = np.empty((n_boot, 2))
    for i in range(n_boot):
        idx = np.concatenate([g[rng.integers(0, len(g), len(g))] for g in groups])
        a, b, _ = fit_scaling(L[idx], V[idx], n[idx])
        out[i] = a, b
    lo, hi = np.percentile(out, [2.5, 97.5], axis=0)
    return (float(lo[0]), float(hi[0])), (float(lo[1]), float(hi[1]))


def mst_tree(points: np.ndarray):
    """Minimum spanning tree of a point cloud (the wiring-optimal tree without
    a root constraint). Returns vertices and edges. Used only for self-tests."""
    from scipy.spatial import Delaunay
    tri = Delaunay(points)
    s = tri.simplices
    pairs = np.vstack([s[:, [i, j]] for i in range(4) for j in range(i + 1, 4)])
    pairs = np.unique(np.sort(pairs, axis=1), axis=0)
    w = np.linalg.norm(points[pairs[:, 0]] - points[pairs[:, 1]], axis=1)
    m = coo_matrix((w, (pairs[:, 0], pairs[:, 1])), shape=(len(points),) * 2)
    t = minimum_spanning_tree(m).tocoo()
    return points, np.column_stack([t.row, t.col]).astype(np.int64)
