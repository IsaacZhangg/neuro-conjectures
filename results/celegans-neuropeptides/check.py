"""Do neuropeptide links improve prediction of C. elegans signal propagation?

Pre-registered test (see CLAIM.md and METHODS.md). Compares a parameter-free
path predictor, S(A) = A + A.A, on the wired connectome and on the wired plus
neuropeptide connectome, scored by AUROC against the wild-type responses of
the Randi et al. 2023 signal-propagation atlas.

    python check.py           full run, writes check-output.json, sets status
    python check.py --quick   smoke test (under a minute), no status

Exit code: full run 0 if the primary analysis is Supported, 1 otherwise.
Quick mode 0 if the pipeline ran.
"""
from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SEED = 20261125
Q_TH = 0.05
RANGES = ["short_range_model", "mid_range_model", "long_range_model"]
PRIMARY = "PRIMARY_NETWORK_PLACEHOLDER"
N_BOOT, N_NULL = 2000, 1000
N_BOOT_QUICK, N_NULL_QUICK = 50, 20
SWAPS_PER_EDGE, MAX_ATTEMPTS_PER_EDGE = 10, 100
PREFIXES = ["DA", "DB", "AS", "VA", "VB", "VC", "VD", "DD"]


# ---------------------------------------------------------------- loading

def standardise(name: str) -> str:
    name = name.strip()
    for pf in PREFIXES:
        name = name.replace(pf + "0", pf)
    return {"AWCL": "AWCOFF", "AWCR": "AWCON"}.get(name, name)


def load_peptide(path: Path, ids: list[str]) -> tuple[np.ndarray, list[str]]:
    """Binary peptide network on the atlas order, indexed [post, pre]."""
    index = {n: k for k, n in enumerate(ids)}
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    header = [standardise(x) for x in rows[0][1:]]
    A = np.zeros((len(ids), len(ids)), dtype=bool)
    unmatched = sorted({n for n in header if n not in index})
    for r in rows[1:]:
        pre = standardise(r[0])
        if pre not in index:
            unmatched.append(pre)
            continue
        j = index[pre]
        for post, v in zip(header, r[1:]):
            if post in index and v.strip() and float(v) > 0:
                A[index[post], j] = True
    np.fill_diagonal(A, False)
    return A, sorted(set(unmatched))


def load_real(data_dir: Path) -> dict:
    import h5py

    ids = list(np.loadtxt(data_dir / "neuron_ids.txt", dtype=str)[:, 1])
    out = {"ids": ids}
    with h5py.File(data_dir / "funatlas.h5", "r") as f:
        h5ids = [n.decode() for n in f["neuron_ids"][:]]
        mismatched = [(a, b) for a, b in zip(h5ids, ids) if a != b]
        out["atlas_id_mismatches"] = mismatched
        for strain in ["wt", "unc31"]:
            out[strain] = {k: f[strain][k][:] for k in ["q", "q_eq", "occ1"]}
    with h5py.File(data_dir / "aconnectome_default.h5", "r") as f:
        chem, gap = f["chem"][:], f["gap"][:]
    W = (chem > 0) | (gap > 0) | (gap.T > 0)
    np.fill_diagonal(W, False)
    out["A_w"] = W
    out["peptide"] = {}
    for rng_name in RANGES:
        A, unmatched = load_peptide(
            data_dir / f"01022024_neuropeptide_connectome_{rng_name}.csv", ids)
        out["peptide"][rng_name] = {"A": A, "unmatched": unmatched}
    return out


def load_synthetic(seed: int) -> dict:
    """Small fake dataset for the CI smoke test when the data are unavailable."""
    rng = np.random.default_rng(seed)
    n = 40
    W = rng.random((n, n)) < 0.05
    P = rng.random((n, n)) < 0.15
    np.fill_diagonal(W, False)
    np.fill_diagonal(P, False)
    drive = W.astype(float) + 0.5 * P + 0.2 * rng.random((n, n))
    q = np.where(drive > 0.6, 0.01, 0.5)
    q_eq = np.where(drive < 0.3, 0.01, 0.5)
    occ = (rng.random((n, n)) < 0.6).astype(int)
    strain = {"q": q, "q_eq": q_eq, "occ1": occ}
    pep = {r: {"A": P.copy(), "unmatched": []} for r in RANGES}
    return {"ids": [f"N{k}" for k in range(n)], "wt": strain, "unc31": strain,
            "A_w": W, "peptide": pep, "atlas_id_mismatches": []}


# ---------------------------------------------------------------- statistics

def path_score(A: np.ndarray, two_step: bool = True) -> np.ndarray:
    A = A.astype(np.int64)
    return A + A @ A if two_step else A


def auroc(score: np.ndarray, y: np.ndarray) -> float:
    """Mann-Whitney AUROC with mid-ranks for ties."""
    n1 = int(y.sum())
    n0 = y.size - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    s = score[order]
    ranks = np.empty(score.size, dtype=float)
    # mid-ranks
    edges = np.flatnonzero(np.diff(s)) + 1
    starts = np.concatenate(([0], edges))
    ends = np.concatenate((edges, [s.size]))
    mid = (starts + ends + 1) / 2.0
    ranks[order] = np.repeat(mid, ends - starts)
    return float((ranks[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def gains(sw: np.ndarray, su: np.ndarray, y: np.ndarray, none: np.ndarray) -> dict:
    g_all = auroc(su, y) - auroc(sw, y)
    g_none = auroc(su[none], y[none]) - 0.5
    path = ~none
    g_path = auroc(su[path], y[path]) - auroc(sw[path], y[path])
    return {"G_all": g_all, "G_none": g_none, "G_path": g_path}


def rewire_py(A: np.ndarray, seed: int) -> np.ndarray:
    """Degree-preserving rewiring (row and column sums kept, no diagonal)."""
    rng = np.random.default_rng(seed)
    M = A.copy()
    r, c = np.nonzero(M)
    e = r.size
    target, attempts, accepted = SWAPS_PER_EDGE * e, 0, 0
    while accepted < target and attempts < MAX_ATTEMPTS_PER_EDGE * e:
        attempts += 1
        x, z = rng.integers(e, size=2)
        r1, c1, r2, c2 = r[x], c[x], r[z], c[z]
        if r1 == r2 or c1 == c2 or r1 == c2 or r2 == c1 or M[r1, c2] or M[r2, c1]:
            continue
        M[r1, c1] = M[r2, c2] = False
        M[r1, c2] = M[r2, c1] = True
        c[x], c[z] = c2, c1
        accepted += 1
    return M


try:
    import numba

    @numba.njit(cache=False)
    def _rewire_nb(M, r, c, seed, target, max_attempts):
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
            M[r1, c1] = False
            M[r2, c2] = False
            M[r1, c2] = True
            M[r2, c1] = True
            c[x] = c2
            c[z] = c1
            accepted += 1
        return accepted

    def rewire(A: np.ndarray, seed: int) -> np.ndarray:
        M = A.copy()
        r, c = np.nonzero(M)
        r, c = r.astype(np.int64), c.astype(np.int64)
        e = r.size
        _rewire_nb(M, r, c, seed % (2**32 - 1), SWAPS_PER_EDGE * e, MAX_ATTEMPTS_PER_EDGE * e)
        return M

    REWIRER = f"numba {numba.__version__}"
except ImportError:  # CI smoke test path
    rewire = rewire_py
    REWIRER = "pure python"


# ---------------------------------------------------------------- analysis

def pair_set(strain: dict, label_mode: str):
    q, q_eq, occ = strain["q"], strain["q_eq"], strain["occ1"]
    n = q.shape[0]
    measured = (occ >= 1) & np.isfinite(q) & ~np.eye(n, dtype=bool)
    pos = measured & (q < Q_TH)
    if label_mode == "standard":
        keep = measured
    elif label_mode == "q_eq_negatives":
        keep = pos | (measured & ~pos & np.isfinite(q_eq) & (q_eq < Q_TH))
    else:
        raise ValueError(label_mode)
    ii, jj = np.nonzero(keep)
    return ii, jj, pos[ii, jj]


def _r(x, nd: int = 5):
    x = float(x)
    return round(x, nd) if np.isfinite(x) else None


def _passes(entry: dict) -> bool:
    lo, p = entry["ci95"][0], entry.get("null_p")
    return lo is not None and p is not None and lo > 0 and p < 0.05


def decide(st: dict) -> dict:
    a = _passes(st["G_all"])
    b = _passes(st["G_none"])
    gn, gp = st["G_none"]["value"], st["G_path"]["value"]
    c = gn is not None and gp is not None and gn > gp
    if a and b and c:
        status = "Supported"
    elif (a and b) or (a != b):
        status = "Partial"
    else:
        status = "Refuted"
    return {"a": bool(a), "b": bool(b), "c": bool(c), "status": status}


def run(data: dict, n_boot: int, n_null: int, variants: list[dict], seed: int) -> dict:
    rng = np.random.default_rng(seed)
    A_w = data["A_w"]
    results = {}
    # set up each variant's pairs and observed values
    prepared = {}
    for v in variants:
        ii, jj, y = pair_set(data["wt"], v["labels"])
        two = v["two_step"]
        sw_full = path_score(A_w, two)
        A_p = data["peptide"][v["network"]]["A"]
        su_full = path_score(A_w | A_p, two)
        sw, su = sw_full[ii, jj], su_full[ii, jj]
        none = sw == 0
        obs = gains(sw, su, y, none)
        groups = {}
        for k, j in enumerate(jj):
            groups.setdefault(int(j), []).append(k)
        group_idx = [np.array(g) for g in groups.values()]
        prepared[v["name"]] = dict(ii=ii, jj=jj, y=y, sw=sw, su=su, none=none, obs=obs,
                                   groups=group_idx, v=v)
    # bootstrap
    for name, p in prepared.items():
        bs = {k: [] for k in p["obs"]}
        g = p["groups"]
        for _ in range(n_boot):
            draw = rng.integers(len(g), size=len(g))
            idx = np.concatenate([g[d] for d in draw])
            gb = gains(p["sw"][idx], p["su"][idx], p["y"][idx], p["none"][idx])
            for k in bs:
                bs[k].append(gb[k])
        p["boot"] = {k: np.array(x) for k, x in bs.items()}
    # null: rewire each network once per null draw, reuse across its variants
    null_seeds = rng.integers(2**31 - 1, size=n_null)
    for net in sorted({p["v"]["network"] for p in prepared.values()}):
        A_p = data["peptide"][net]["A"]
        users = [p for p in prepared.values() if p["v"]["network"] == net]
        for p in users:
            p["null"] = {"G_all": [], "G_none": []}
        changed = []
        for s in null_seeds:
            R = rewire(A_p, int(s))
            changed.append(float((R & ~A_p).sum() / max(1, A_p.sum())))
            for p in users:
                two = p["v"]["two_step"]
                su = path_score(A_w | R, two)[p["ii"], p["jj"]]
                gn = gains(p["sw"], su, p["y"], p["none"])
                p["null"]["G_all"].append(gn["G_all"])
                p["null"]["G_none"].append(gn["G_none"])
        for p in users:
            p["null_fraction_edges_moved"] = float(np.mean(changed)) if changed else None
    for name, p in prepared.items():
        st = {}
        for k in ["G_all", "G_none", "G_path"]:
            b = p["boot"][k]
            b = b[np.isfinite(b)]
            entry = {"value": _r(p["obs"][k]),
                     "ci95": [_r(np.percentile(b, 2.5)), _r(np.percentile(b, 97.5))] if b.size else [None, None],
                     "bootstrap_finite": int(b.size)}
            if k in p["null"]:
                nul = np.array(p["null"][k])
                nul = nul[np.isfinite(nul)]
                ok = nul.size and np.isfinite(p["obs"][k])
                entry["null_p"] = _r((1 + np.sum(nul >= p["obs"][k])) / (1 + nul.size)) if ok else None
                entry["null_mean"] = _r(nul.mean()) if nul.size else None
                entry["null_95th"] = _r(np.percentile(nul, 95)) if nul.size else None
            st[k] = entry
        y, none = p["y"], p["none"]
        st["auroc_wired_all"] = _r(auroc(p["sw"], y))
        st["auroc_combined_all"] = _r(auroc(p["su"], y))
        st["counts"] = {
            "pairs": int(y.size), "responding": int(y.sum()),
            "stimulated_neurons": len(p["groups"]),
            "pairs_no_wired_path": int(none.sum()),
            "responding_no_wired_path": int(y[none].sum()),
        }
        st["null_mean_fraction_peptide_edges_moved"] = p["null_fraction_edges_moved"]
        st["decision"] = decide(st)
        st["spec"] = p["v"]
        results[name] = st
    return results


def exploratory_unc31(data: dict, net: str) -> dict:
    A_w, A_p = data["A_w"], data["peptide"][net]["A"]
    none = path_score(A_w) == 0
    out = {}
    for strain in ["wt", "unc31"]:
        ii, jj, y = pair_set(data[strain], "standard")
        sel = none[ii, jj]
        pep = A_p[ii, jj]
        row = {}
        for lab, m in [("peptide_link", sel & pep), ("no_peptide_link", sel & ~pep)]:
            row[lab] = {"pairs": int(m.sum()), "responding": int(y[m].sum()),
                        "rate": round(float(y[m].mean()), 4) if m.any() else None}
        out[strain] = row
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--data", default=str(HERE / "data"))
    args = ap.parse_args()
    t0 = time.time()

    source = "real"
    try:
        sys.path.insert(0, str(HERE))
        from fetch_data import fetch
        checksums = fetch(Path(args.data))
        data = load_real(Path(args.data))
    except Exception as exc:  # noqa: BLE001
        if not args.quick:
            print(f"could not load data: {exc}")
            return 1
        print(f"data unavailable ({type(exc).__name__}: {exc}); running synthetic smoke test")
        source, checksums, data = "synthetic", {}, load_synthetic(args.seed)

    if PRIMARY not in RANGES:
        raise SystemExit("primary network not set in check.py")
    others = [r for r in RANGES if r != PRIMARY]
    variants = [{"name": "primary", "network": PRIMARY, "labels": "standard", "two_step": True}]
    if not args.quick:
        variants += [{"name": f"R1_{r}", "network": r, "labels": "standard", "two_step": True} for r in others]
        variants += [{"name": "R2_q_eq_negatives", "network": PRIMARY, "labels": "q_eq_negatives", "two_step": True},
                     {"name": "R3_direct_only", "network": PRIMARY, "labels": "standard", "two_step": False}]
    n_boot, n_null = (N_BOOT_QUICK, N_NULL_QUICK) if args.quick else (N_BOOT, N_NULL)
    res = run(data, n_boot, n_null, variants, args.seed)

    for name, st in res.items():
        d = st["decision"]
        print(f"{name:22s} G_all={st['G_all']['value']} {st['G_all']['ci95']} p={st['G_all']['null_p']}"
              f"  G_none={st['G_none']['value']} {st['G_none']['ci95']} p={st['G_none']['null_p']}"
              f"  G_path={st['G_path']['value']}  -> {d['status']}")

    if args.quick:
        print(f"quick smoke test on {source} data: pipeline ran; status not set by quick mode")
        return 0

    primary = res["primary"]["decision"]["status"]
    poles = {"Supported": "Refuted", "Refuted": "Supported"}
    flipped = [n for n, st in res.items() if n != "primary"
               and poles.get(primary) == st["decision"]["status"]]
    status = "Contested" if flipped else primary
    out = {
        "result": "celegans-neuropeptides",
        "mode": "full",
        "seed": args.seed,
        "primary_network": PRIMARY,
        "q_threshold": Q_TH,
        "n_bootstrap": n_boot,
        "n_null": n_null,
        "status": status,
        "variants_opposite_to_primary": flipped,
        "analyses": res,
        "exploratory_unc31_no_wired_path": exploratory_unc31(data, PRIMARY),
        "unmatched_peptide_neurons": {r: data["peptide"][r]["unmatched"] for r in RANGES},
        "peptide_links": {r: int(data["peptide"][r]["A"].sum()) for r in RANGES},
        "wired_links": int(data["A_w"].sum()),
        "atlas_neurons": len(data["ids"]),
        "atlas_id_mismatches": data["atlas_id_mismatches"],
        "data_sha256": checksums,
        "rewiring": REWIRER,
        "runtime_seconds": round(time.time() - t0, 1),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
    }
    with open(HERE / "check-output.json", "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
    print(f"status: {status}")
    return 0 if status == "Supported" else 1


if __name__ == "__main__":
    raise SystemExit(main())
