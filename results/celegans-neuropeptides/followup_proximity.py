"""Follow-up: does physical proximity explain the pair-specific peptide gain?

Plan: METHODS-followup-proximity.md (frozen before this was run).

    python followup_proximity.py           full run, writes followup-proximity-output.json
    python followup_proximity.py --quick   smoke test (few resamples)

Does not set the result's status. Exit 0 if it ran.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check  # noqa: E402
from fetch_data import fetch  # noqa: E402

N_BOOT, N_NULL = 2000, 1000
ATTEMPTS_PER_EDGE_F2 = 1000


def load_proximity(data_dir: Path, ids: list[str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    index = {n: k for k, n in enumerate(ids)}
    lines = (data_dir / "anatlas_neuron_positions.txt").read_text().splitlines()
    names = [check.standardise(x) for x in lines[0].lstrip("#").split()]
    coords = np.array([[float(v) for v in ln.split()] for ln in lines[1:1 + len(names)]])
    pos = np.full((len(ids), 3), np.nan)
    missing = []
    for name, xyz in zip(names, coords):
        if name in index:
            pos[index[name]] = xyz
    missing = [n for n in ids if np.isnan(pos[index[n], 0])]
    D = np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(-1))
    gang = json.loads((data_dir / "aconnectome_ids_ganglia.json").read_text())
    member = {}
    for g, neurons in gang.items():
        if g in ("head", "pharynx"):
            continue
        for n in neurons:
            member.setdefault(check.standardise(n), set()).add(g)
    n = len(ids)
    G = np.zeros((n, n), dtype=bool)
    for a in range(n):
        ma = member.get(ids[a], set())
        if not ma:
            continue
        for b in range(n):
            if ma & member.get(ids[b], set()):
                G[a, b] = True
    return D, G, missing


def proximity_classes(D: np.ndarray, G: np.ndarray) -> tuple[np.ndarray, list[float]]:
    n = D.shape[0]
    off = ~np.eye(n, dtype=bool)
    ok = np.isfinite(D) & off
    edges = list(np.quantile(D[ok], [0.2, 0.4, 0.6, 0.8]))
    q = np.digitize(D, edges)
    C = np.where(ok, q * 2 + G.astype(int), -1).astype(np.int64)
    return C, [float(e) for e in edges]


try:
    import numba

    @numba.njit(cache=False)
    def _rewire_class_nb(M, r, c, C, seed, target, max_attempts):
        np.random.seed(seed)
        e = r.size
        attempts = 0
        accepted = 0
        while accepted < target and attempts < max_attempts:
            attempts += 1
            x = np.random.randint(e)
            z = np.random.randint(e)
            r1, c1, r2, c2 = r[x], c[x], r[z], c[z]
            if r1 == r2 or c1 == c2 or r1 == c2 or r2 == c1:
                continue
            if M[r1, c2] or M[r2, c1]:
                continue
            if C[r1, c2] != C[r1, c1] or C[r2, c1] != C[r2, c2]:
                continue
            M[r1, c1] = False
            M[r2, c2] = False
            M[r1, c2] = True
            M[r2, c1] = True
            c[x] = c2
            c[z] = c1
            accepted += 1
        return accepted

    def rewire_class(A, C, seed):
        M = A.copy()
        r, c = np.nonzero(M)
        r, c = r.astype(np.int64), c.astype(np.int64)
        e = r.size
        acc = _rewire_class_nb(M, r, c, C, seed % (2**32 - 1),
                               check.SWAPS_PER_EDGE * e, ATTEMPTS_PER_EDGE_F2 * e)
        return M, acc
except ImportError:
    def rewire_class(A, C, seed):
        rng = np.random.default_rng(seed)
        M = A.copy()
        r, c = np.nonzero(M)
        e = r.size
        acc = att = 0
        while acc < check.SWAPS_PER_EDGE * e and att < ATTEMPTS_PER_EDGE_F2 * e:
            att += 1
            x, z = rng.integers(e, size=2)
            r1, c1, r2, c2 = r[x], c[x], r[z], c[z]
            if r1 == r2 or c1 == c2 or r1 == c2 or r2 == c1 or M[r1, c2] or M[r2, c1]:
                continue
            if C[r1, c2] != C[r1, c1] or C[r2, c1] != C[r2, c2]:
                continue
            M[r1, c1] = M[r2, c2] = False
            M[r1, c2] = M[r2, c1] = True
            c[x], c[z] = c2, c1
            acc += 1
        return M, acc


def summarise(real: dict, boot: dict, null: dict) -> dict:
    out = {}
    for k in ["G_all", "G_none"]:
        nul = np.array(null[k])
        b = np.array(boot[k])
        b = b[np.isfinite(b)]
        m = float(nul.mean())
        out[k] = {
            "value": check._r(real[k]),
            "ci95": [check._r(np.percentile(b, 2.5)), check._r(np.percentile(b, 97.5))],
            "null_mean": check._r(m),
            "null_95th": check._r(np.percentile(nul, 95)),
            "null_p": check._r((1 + np.sum(nul >= real[k])) / (1 + nul.size)),
            "excess": check._r(real[k] - m),
            "excess_ci95": [check._r(np.percentile(b, 2.5) - m), check._r(np.percentile(b, 97.5) - m)],
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=check.SEED)
    ap.add_argument("--data", default=str(HERE / "data"))
    args = ap.parse_args()
    t0 = time.time()
    n_boot, n_null = (50, 20) if args.quick else (N_BOOT, N_NULL)
    fetch(Path(args.data))
    data = check.load_real(Path(args.data))
    ids = data["ids"]
    D, Gang, missing = load_proximity(Path(args.data), ids)
    C, q_edges = proximity_classes(D, Gang)
    A_w, A_p = data["A_w"], data["peptide"][check.PRIMARY]["A"]

    ii, jj, y, _ = check.pair_set(data["wt"], "standard")
    sw = check.path_score(A_w)[ii, jj]
    su = check.path_score(A_w | A_p)[ii, jj]
    none = sw == 0
    strata = C[ii, jj]
    groups = {}
    for k, j in enumerate(jj):
        groups.setdefault(int(j), []).append(k)
    groups = [np.array(g) for g in groups.values()]
    rng = np.random.default_rng(args.seed)

    def boot(st):
        res = {"G_all": [], "G_none": []}
        for _ in range(n_boot):
            idx = np.concatenate([groups[d] for d in rng.integers(len(groups), size=len(groups))])
            g = check.gains(sw[idx], su[idx], y[idx], none[idx], None if st is None else st[idx])
            res["G_all"].append(g["G_all"])
            res["G_none"].append(g["G_none"])
        return res

    # F1: stratified by proximity class, ordinary degree-preserving null
    real1 = check.gains(sw, su, y, none, strata)
    boot1 = boot(strata)
    null1 = {"G_all": [], "G_none": []}
    seeds1 = rng.integers(2**31 - 1, size=n_null)
    for s in seeds1:
        R = check.rewire(A_p, int(s))
        g = check.gains(sw, check.path_score(A_w | R)[ii, jj], y, none, strata)
        null1["G_all"].append(g["G_all"])
        null1["G_none"].append(g["G_none"])

    # F2: unstratified, proximity-class-preserving null
    real2 = check.gains(sw, su, y, none, None)
    boot2 = boot(None)
    null2 = {"G_all": [], "G_none": []}
    moved, accepted = [], []
    seeds2 = rng.integers(2**31 - 1, size=n_null)
    for s in seeds2:
        R, acc = rewire_class(A_p, C, int(s))
        moved.append(float((R & ~A_p).sum() / A_p.sum()))
        accepted.append(int(acc))
        g = check.gains(sw, check.path_score(A_w | R)[ii, jj], y, none, None)
        null2["G_all"].append(g["G_all"])
        null2["G_none"].append(g["G_none"])
        assert (R.sum(0) == A_p.sum(0)).all() and (R.sum(1) == A_p.sum(1)).all()

    f1 = summarise(real1, boot1, null1)
    f2 = summarise(real2, boot2, null2)
    s1 = f1["G_none"]["null_p"] < 0.05
    s2 = f2["G_none"]["null_p"] < 0.05
    verdict = "Survives proximity" if s1 and s2 else ("Mixed" if s1 != s2 else "Explained by proximity")

    # exploratory
    nd = -np.where(np.isfinite(D), D, np.nanmax(D))[ii, jj]
    has_pos = np.isfinite(D[ii, jj])
    off = ~np.eye(len(ids), dtype=bool)
    qbin = np.digitize(D, q_edges)
    share = {}
    for b in range(5):
        m = off & np.isfinite(D) & (qbin == b)
        share[f"quintile_{b + 1}"] = check._r(A_p[m].mean(), 4)
    out = {
        "result": "celegans-neuropeptides follow-up (proximity)",
        "mode": "quick" if args.quick else "full",
        "plan": "METHODS-followup-proximity.md",
        "seed": args.seed, "n_bootstrap": n_boot, "n_null": n_null,
        "network": check.PRIMARY,
        "F1_proximity_stratified": f1,
        "F2_proximity_preserving_shuffle": {
            **f2,
            "mean_fraction_links_moved": check._r(np.mean(moved), 4),
            "mean_accepted_swaps": float(np.mean(accepted)),
            "target_swaps": int(check.SWAPS_PER_EDGE * A_p.sum()),
        },
        "survives_F1": bool(s1), "survives_F2": bool(s2), "verdict": verdict,
        "exploratory": {
            "auroc_nearness_all_pairs_with_positions": check._r(check.auroc(nd[has_pos], y[has_pos])),
            "auroc_nearness_no_wired_path": check._r(check.auroc(nd[has_pos & none], y[has_pos & none])),
            "peptide_link_share_by_distance_quintile": share,
            "auroc_same_ganglion_no_wired_path": check._r(check.auroc(Gang[ii, jj][none].astype(int), y[none])),
        },
        "neurons_without_position": missing,
        "pairs_in_class_minus1": int((strata == -1).sum()),
        "distance_quintile_edges": q_edges,
        "runtime_seconds": round(time.time() - t0, 1),
        "python": sys.version.split()[0], "numpy": np.__version__,
        "platform": platform.platform(),
    }
    for name, f in [("F1", f1), ("F2", f2)]:
        print(name, "G_none", f["G_none"], "G_all", f["G_all"]["value"], f["G_all"]["null_p"])
    print("verdict:", verdict)
    if not args.quick:
        with open(HERE / "followup-proximity-output.json", "w") as fh:
            json.dump(out, fh, indent=2)
            fh.write("\n")
    else:
        print(json.dumps(out["exploratory"]), out["F2_proximity_preserving_shuffle"]["mean_fraction_links_moved"], missing, out["pairs_in_class_minus1"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
