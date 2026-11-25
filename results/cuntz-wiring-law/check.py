"""Cuntz wiring law (L ∝ V^(1/3) bp^(2/3)) on FlyWire and MICrONS skeletons.

    python check.py --quick   # estimator benchmark only, no network (< 1 min)
    python check.py           # full run: downloads skeletons, applies the
                              # frozen decision rule in CLAIM.md

Exit 0 when the result is not Refuted and the benchmark passes, 1 otherwise.
The full run writes check-output.json, neuron-stats.csv and manifest.json.
Downloaded files are cached in ./data/ (not committed).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import platform
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from wiring import fit_scaling, mst_tree, tree_stats

HERE = Path(__file__).resolve().parent
SEED = 20261125
W = (2 / 3 - 0.1, 2 / 3 + 0.1)
TAUS = (0.0, 1.0, 5.0)
N_BOOT = 2000
MIN_L, MIN_BP = 100.0, 10

BOSS = "https://bossdb-open-data.s3.amazonaws.com/iarpa_microns/minnie/minnie65"
MICRONS_SKEL = BOSS + "/skeletons/v661/skeletons/"
MICRONS_TYPES = BOSS + "/cell_types/allen_visp_column_soma_coarse_types_v1.csv"
MICRONS_NUCLEI = BOSS + "/nucleus_detection/nucleus_detection_v0.csv"
FLYWIRE_CV = "precomputed://https://storage.googleapis.com/flywire_v141_m783"
FLYWIRE_ANN_COMMIT = "a83b2776d60d5764cef36b927f5f9679c16c47a2"
FLYWIRE_ANN = ("https://raw.githubusercontent.com/flyconnectome/flywire_annotations/"
               f"{FLYWIRE_ANN_COMMIT}/supplemental_files/Supplemental_file1_neuron_annotations.tsv")
FLY_SUPER = {"central", "optic", "visual_projection", "visual_centrifugal", "endocrine"}
FLY_PER_TYPE = 15

DATA = HERE / "data"


def sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def fetch(url: str, cache: Path, retries: int = 4) -> bytes:
    if cache.exists():
        return cache.read_bytes()
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                b = r.read()
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_bytes(b)
            return b
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(2 ** i)


# ---------------------------------------------------------------- statistics

def beta_fit(L, V, n):
    y = np.log(L / np.cbrt(V))
    return float(np.polyfit(np.log(n), y, 1)[0])


def bootstrap(L, V, n, strata, rng, n_boot=N_BOOT):
    L, V, n, strata = map(np.asarray, (L, V, n, strata))
    groups = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    betas, ab = np.empty(n_boot), np.empty((n_boot, 2))
    for i in range(n_boot):
        idx = np.concatenate([g[rng.integers(0, len(g), len(g))] for g in groups])
        betas[i] = beta_fit(L[idx], V[idx], n[idx])
        a, b, _ = fit_scaling(L[idx], V[idx], n[idx])
        ab[i] = a, b
    ci = lambda x: [float(v) for v in np.percentile(x, [2.5, 97.5])]
    a, b, _ = fit_scaling(L, V, n)
    return {
        "beta": beta_fit(L, V, n), "beta_ci95": ci(betas),
        "free_a": a, "free_a_ci95": ci(ab[:, 0]),
        "free_b": b, "free_b_ci95": ci(ab[:, 1]),
        "n_neurons": int(len(L)),
    }


def classify(ci):
    lo, hi = ci
    if W[0] <= lo and hi <= W[1]:
        return "inside"
    if hi < W[0] or lo > W[1]:
        return "outside"
    return "overlaps"


def dataset_verdict(primary, variants):
    p = classify(primary["beta_ci95"])
    v = [classify(x["beta_ci95"]) for x in variants]
    if p == "inside":
        return "Contested" if any(x == "outside" for x in v) else "Supported"
    if p == "outside":
        if all(x == "outside" for x in v):
            return "Refuted"
        return "Contested" if any(x == "inside" for x in v) else "Inconclusive"
    return "Inconclusive"


def overall_status(verdicts):
    vs = list(verdicts)
    if all(v == "Refuted" for v in vs):
        return "Refuted"
    if any(v == "Contested" for v in vs):
        return "Contested"
    return "Partial"


# ---------------------------------------------------------------- benchmark

def benchmark():
    rng = np.random.default_rng(SEED)
    rows = []
    for _ in range(150):
        npts = int(np.exp(rng.uniform(np.log(50), np.log(3000))))
        side = np.exp(rng.uniform(np.log(50), np.log(500)))
        shape = np.array([1, rng.uniform(0.3, 1), rng.uniform(0.3, 1)]) * side
        s = tree_stats(*mst_tree(rng.uniform(0, 1, (npts, 3)) * shape), 0.0)
        rows.append((s.L, s.V_hull, s.n_branch))
    r = np.array(rows)
    beta = beta_fit(r[:, 0], r[:, 1], r[:, 2])
    # pruning self-test on a hand-built tree: a 10 um trunk, a 4 um branch and a 0.5 um twig
    xyz = np.array([[0, 0, 0], [5, 0, 0], [10, 0, 0], [5, 0.5, 0], [5, 0, 4]], float)
    edges = np.array([[0, 1], [1, 2], [1, 3], [1, 4]])
    s0, s1 = tree_stats(xyz, edges, 0.0), tree_stats(xyz, edges, 1.0)
    prune_ok = (s0.n_branch == 1 and s0.n_tips == 4 and abs(s0.L - 14.5) < 1e-9
                and s1.n_branch == 1 and s1.n_tips == 3 and abs(s1.L - 14.0) < 1e-9)
    return {"mst_beta": beta, "mst_trees": 150, "window": list(W),
            "mst_inside_window": bool(W[0] <= beta <= W[1]), "prune_selftest": bool(prune_ok)}


# ---------------------------------------------------------------- MICrONS

def parse_swc(b: bytes):
    a = np.loadtxt(io.StringIO(b.decode()), comments="#", ndmin=2)
    ids, typ, xyz, par = a[:, 0].astype(int), a[:, 1].astype(int), a[:, 2:5], a[:, 6].astype(int)
    keep = np.isin(typ, [1, 3, 4])
    index = {i: k for k, i in enumerate(ids[keep])}
    e = [(index[i], index[p]) for i, p in zip(ids[keep], par[keep]) if p in index]
    return xyz[keep], np.array(e, dtype=np.int64).reshape(-1, 2)


def microns(manifest, quick_n=None):
    types_b = fetch(MICRONS_TYPES, DATA / "microns" / "coarse_types.csv")
    manifest["microns_cell_types_sha256"] = sha256(types_b)
    # The column table's first field is an annotation id, not a nucleus id
    # (deviation 1 in METHODS.md): map each row to a nucleus by exact soma
    # position in nucleus_detection_v0.csv.
    nd_b = fetch(MICRONS_NUCLEI, DATA / "microns" / "nucleus_detection_v0.csv")
    manifest["microns_nucleus_table_sha256"] = sha256(nd_b)
    pos2nuc = {}
    for row in csv.reader(io.StringIO(nd_b.decode())):
        pos2nuc.setdefault((row[4], row[5], row[6]), []).append(int(row[0]))
    cls, unmatched = {}, 0
    for row in csv.reader(io.StringIO(types_b.decode())):
        if row[1] == "t" and row[3] in ("excitatory", "inhibitory"):
            hit = pos2nuc.get((row[4], row[5], row[6]), [])
            if len(hit) == 1:
                cls[hit[0]] = row[3]
            else:
                unmatched += 1
    keys = list_s3(MICRONS_SKEL)
    by_nuc = {}
    for k in keys:
        stem = k.rsplit("/", 1)[-1][:-4]
        root, nuc = stem.split("_")
        by_nuc.setdefault(int(nuc), []).append(k)
    todo, excl = [], {"no_unique_nucleus_match": unmatched, "no_skeleton": 0, "ambiguous_skeleton": 0}
    for nuc in sorted(cls):
        ks = by_nuc.get(nuc, [])
        if len(ks) == 0:
            excl["no_skeleton"] += 1
        elif len(ks) > 1:
            excl["ambiguous_skeleton"] += 1
        else:
            todo.append((nuc, ks[0]))
    if quick_n:
        todo = todo[:quick_n]

    def work(item):
        nuc, key = item
        b = fetch("https://bossdb-open-data.s3.amazonaws.com/" + key, DATA / "microns" / key.rsplit("/", 1)[-1])
        xyz, e = parse_swc(b)
        return nuc, key, sha256(b), {t: tree_stats(xyz, e, t) for t in TAUS}

    with ThreadPoolExecutor(16) as ex:
        out = list(ex.map(work, todo))
    rows, files = [], {}
    for nuc, key, h, st in out:
        files[key.rsplit("/", 1)[-1]] = h
        rows.append({"dataset": "MICrONS", "id": key.rsplit("/", 1)[-1][:-4], "stratum": cls[nuc], **flatten(st)})
    manifest["microns_skeleton_sha256"] = files
    manifest["microns_skeletons_listed"] = len(keys)
    return rows, excl


def list_s3(prefix_url):
    base, prefix = prefix_url.split(".com/", 1)
    cache = DATA / "microns" / "listing.txt"
    if cache.exists():
        return cache.read_text().split()
    keys, token = [], None
    import re, urllib.parse
    while True:
        url = f"{base}.com/?list-type=2&prefix={urllib.parse.quote(prefix)}"
        if token:
            url += "&continuation-token=" + urllib.parse.quote(token)
        with urllib.request.urlopen(url, timeout=120) as r:
            body = r.read().decode()
        keys += re.findall(r"<Key>([^<]+\.swc)</Key>", body)
        m = re.search(r"<NextContinuationToken>([^<]+)<", body)
        if not m:
            break
        token = m.group(1)
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text("\n".join(keys))
    return keys


# ---------------------------------------------------------------- FlyWire

def flywire(manifest, quick_n=None):
    import pandas as pd
    from cloudvolume import CloudVolume

    ann_b = fetch(FLYWIRE_ANN, DATA / "flywire" / "annotations.tsv")
    manifest["flywire_annotations_sha256"] = sha256(ann_b)
    a = pd.read_csv(io.BytesIO(ann_b), sep="\t", low_memory=False,
                    usecols=["root_id", "super_class", "cell_type", "status"])
    a = a[a.super_class.isin(FLY_SUPER) & a.status.isna() & a.cell_type.notna()]
    counts = a.cell_type.value_counts()
    types = sorted(counts[counts >= FLY_PER_TYPE].index)
    rng = np.random.default_rng(SEED)
    todo = []
    for t in types:
        ids = np.sort(a.loc[a.cell_type == t, "root_id"].astype("int64").to_numpy())
        todo += [(int(i), t) for i in rng.choice(ids, FLY_PER_TYPE, replace=False)]
    manifest["flywire_types"] = len(types)
    if quick_n:
        todo = todo[:quick_n]
    cachef = DATA / "flywire" / "stats.jsonl"
    done = {}
    if cachef.exists():
        for line in cachef.read_text().splitlines():
            d = json.loads(line)
            done[d["id"]] = d

    def work(item):
        rid, t = item
        if str(rid) in done:
            return done[str(rid)]
        for i in range(4):
            try:
                cv = CloudVolume(FLYWIRE_CV, use_https=True, progress=False)
                s = cv.skeleton.get(rid)
                break
            except Exception as exc:
                if i == 3:
                    return {"id": str(rid), "stratum": t, "error": repr(exc)[:200]}
                time.sleep(2 ** i)
        xyz = s.vertices.astype(np.float64) / 1000.0
        e = s.edges.astype(np.int64)
        h = sha256(s.vertices.astype(np.float32).tobytes() + s.edges.astype(np.uint32).tobytes())
        st = {t_: tree_stats(xyz, e, t_) for t_ in TAUS}
        return {"id": str(rid), "stratum": t, "sha256": h, **flatten(st)}

    rows = []
    with ThreadPoolExecutor(8) as ex, cachef.open("a") as fh:
        for d in ex.map(work, todo):
            if str(d["id"]) not in done:
                fh.write(json.dumps(d) + "\n")
                fh.flush()
            rows.append(d)
    excl = {"read_error": sum("error" in d for d in rows)}
    manifest["flywire_skeleton_sha256"] = {d["id"]: d["sha256"] for d in rows if "sha256" in d}
    rows = [{"dataset": "FlyWire", **{k: v for k, v in d.items() if k != "sha256"}} for d in rows if "error" not in d]
    return rows, excl


def flatten(stats):
    out = {}
    for t, s in stats.items():
        k = f"tau{t:g}"
        out.update({f"L_{k}": s.L, f"bp_{k}": s.n_branch, f"tips_{k}": s.n_tips,
                    f"V_{k}": s.V_hull, f"components_{k}": s.n_components})
    return out


# ---------------------------------------------------------------- analysis

def analyse(rows, excl, rng):
    def arr(key):
        return np.array([r[key] for r in rows], float)

    ok = (arr("L_tau1") >= MIN_L) & (arr("bp_tau1") >= MIN_BP) & (arr("V_tau1") > 0)
    # the same neurons enter every variant; a variant drops a neuron only if
    # its own n or V is zero (log undefined)
    excl = dict(excl, failed_size_rule=int((~ok).sum()))
    sel = [r for r, k in zip(rows, ok) if k]
    strata = np.array([r["stratum"] for r in sel])

    def run(tau, nkey):
        L = np.array([r[f"L_tau{tau}"] for r in sel], float)
        V = np.array([r[f"V_tau{tau}"] for r in sel], float)
        n = np.array([r[f"{nkey}_tau{tau}"] for r in sel], float)
        m = (n > 0) & (V > 0)
        return bootstrap(L[m], V[m], n[m], strata[m], rng)

    res = {"primary": run("1", "bp"), "R1_tau0": run("0", "bp"),
           "R2_tau5": run("5", "bp"), "R3_tips": run("1", "tips")}
    res["verdict"] = dataset_verdict(res["primary"], [res["R1_tau0"], res["R2_tau5"], res["R3_tips"]])
    # secondary: per-stratum beta and the constant
    per = {}
    for s in np.unique(strata):
        m = strata == s
        L = np.array([r["L_tau1"] for r in sel])[m]
        V = np.array([r["V_tau1"] for r in sel])[m]
        n = np.array([r["bp_tau1"] for r in sel])[m]
        if m.sum() >= 10 and np.ptp(np.log(n)) > 0:
            per[str(s)] = beta_fit(L, V, n)
    vals = np.array(list(per.values()))
    L = np.array([r["L_tau1"] for r in sel]); V = np.array([r["V_tau1"] for r in sel])
    n = np.array([r["bp_tau1"] for r in sel])
    res["secondary"] = {
        "per_stratum_beta": per,
        "per_stratum_beta_mean": float(vals.mean()) if len(vals) else None,
        "per_stratum_beta_sd": float(vals.std(ddof=1)) if len(vals) > 1 else None,
        "per_stratum_beta_share_inside_window": float(np.mean((vals >= W[0]) & (vals <= W[1]))) if len(vals) else None,
        "constant_median_L_over_V13_bp23": float(np.median(L / (np.cbrt(V) * n ** (2 / 3)))),
        "median_L_um": float(np.median(L)), "median_bp": float(np.median(n)),
        "median_V_um3": float(np.median(V)),
    }
    res["exclusions"] = excl
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="benchmark only, no network")
    ap.add_argument("--smoke", type=int, default=0, help="full pipeline on the first N neurons of each dataset (not a result)")
    args = ap.parse_args()
    t0 = time.time()
    out = {"claim": "L ∝ V^(1/3) · bp^(2/3) (Cuntz, Mathy & Häusser 2012)", "seed": SEED,
           "mode": "quick" if args.quick else ("smoke" if args.smoke else "full"),
           "versions": {"python": platform.python_version(), "numpy": np.__version__}}
    import scipy
    out["versions"]["scipy"] = scipy.__version__
    out["benchmark"] = bench = benchmark()
    print(f"benchmark: MST beta = {bench['mst_beta']:.3f}, window {W[0]:.3f}-{W[1]:.3f}, prune self-test {bench['prune_selftest']}")
    if not (bench["mst_inside_window"] and bench["prune_selftest"]):
        print("benchmark failed: the estimator cannot recognise wiring-optimal trees")
        return 1
    if args.quick:
        print(f"quick mode done in {time.time() - t0:.1f} s (smoke test only; sets no status)")
        return 0
    import cloudvolume
    import pandas
    out["versions"].update({"pandas": pandas.__version__, "cloud-volume": cloudvolume.__version__})
    manifest = {"microns_skeletons": MICRONS_SKEL, "microns_cell_types": MICRONS_TYPES,
                "flywire_skeletons": FLYWIRE_CV + " (skeletons_mip_1)", "flywire_annotations": FLYWIRE_ANN,
                "download_date_utc": time.strftime("%Y-%m-%d", time.gmtime())}
    rng = np.random.default_rng(SEED)
    q = args.smoke or None
    m_rows, m_excl = microns(manifest, q)
    print(f"MICrONS: {len(m_rows)} skeletons measured")
    f_rows, f_excl = flywire(manifest, q)
    print(f"FlyWire: {len(f_rows)} skeletons measured")
    out["datasets"] = {"MICrONS": analyse(m_rows, m_excl, rng), "FlyWire": analyse(f_rows, f_excl, rng)}
    out["status"] = overall_status(d["verdict"] for d in out["datasets"].values())
    out["runtime_s"] = round(time.time() - t0, 1)
    suffix = "" if not args.smoke else "-smoke"
    (HERE / f"check-output{suffix}.json").write_text(json.dumps(out, indent=2) + "\n")
    (HERE / f"manifest{suffix}.json").write_text(json.dumps(manifest, indent=1) + "\n")
    with (HERE / f"neuron-stats{suffix}.csv").open("w", newline="") as fh:
        rows = m_rows + f_rows
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for k, d in out["datasets"].items():
        p = d["primary"]
        print(f"{k}: beta = {p['beta']:.3f} [{p['beta_ci95'][0]:.3f}, {p['beta_ci95'][1]:.3f}], "
              f"n = {p['n_neurons']}, verdict {d['verdict']}")
    print(f"Status: {out['status']}  ({out['runtime_s']} s)")
    return 1 if out["status"] == "Refuted" else 0


if __name__ == "__main__":
    sys.exit(main())
