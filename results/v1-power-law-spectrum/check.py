"""Power-law eigenspectrum of mouse visual cortex responses (alpha ~ 1).

Implements METHODS.md and the decision rule in CLAIM.md.

    python check.py           full run on Allen Brain Observatory data
    python check.py --quick   synthetic end-to-end smoke test (offline, <1 min)

Exit 0 when the overall verdict is "Partial (consistent)" (full run) or when the
synthetic recovery test passes (quick run); exit 1 otherwise. A full run writes
check-output.json next to this file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
S3 = "https://allen-brain-observatory.s3.amazonaws.com/"
P = 118                      # natural images
NR = P - 1                   # usable ranks after centring
AREAS = ["VISp", "VISl", "VISal", "VISrl", "VISpm", "VISam"]
EXCITATORY = ["Cux2", "Emx1", "Slc17a7", "Rorb", "Scnn1a", "Nr5a1", "Rbp4", "Fezf2", "Tlx3"]
B_BAND = (1.00, 1.10)
T_BAND = (0.90, 1.20)
EXACT_DIMS = 3000            # ranks simulated exactly in the forward model
N_SPLITS_MAIN = 10           # trial splits for the point estimate
N_SPLITS_VAR = 4             # trial splits stored for robustness variants
N_BOOT = 500
N_CV = 100
FIT_RANKS = 100

CACHE = Path(os.environ.get("NC_CACHE", Path.home() / ".cache" / "neuro-conjectures" / "v1-power-law-spectrum"))


# --------------------------------------------------------------------------
# cvPCA on Gram matrices
#
# For a pseudo-population X1, X2 (images x neurons, centred over images), the
# cvPCA variance along the k-th principal direction of X1 is
#   lambda_k = u_k^T (X2 X1^T) u_k / (P - 1),  u_k eigenvectors of X1 X1^T.
# All three Gram matrices are sums over neurons, so per-session Grams add up.
# --------------------------------------------------------------------------

def session_grams(resp, frame, rng, n_splits, zscore=True, neuron_mask=None):
    """resp: trials x neurons; frame: image index per trial (0..117, -1 blank).

    Returns array (n_splits, 3, P, P) of [X1X1', X2X1', X2X2'] Gram matrices,
    and the number of neurons used.
    """
    keep = frame >= 0
    resp = resp[keep].astype(np.float64)
    frame = frame[keep]
    if neuron_mask is not None:
        resp = resp[:, neuron_mask]
    if zscore:
        sd = resp.std(0)
        good = sd > 0
        resp = resp[:, good] / sd[good]
    else:
        good = resp.std(0) > 0
        resp = resp[:, good]
    n = resp.shape[1]
    idx = [np.flatnonzero(frame == i) for i in range(P)]
    out = np.zeros((n_splits, 3, P, P))
    for s in range(n_splits):
        X1 = np.empty((P, n))
        X2 = np.empty((P, n))
        for i, ii in enumerate(idx):
            perm = rng.permutation(ii)
            h = len(perm) // 2
            X1[i] = resp[perm[:h]].mean(0)
            X2[i] = resp[perm[h:2 * h]].mean(0)
        X1 -= X1.mean(0)
        X2 -= X2.mean(0)
        out[s, 0] = X1 @ X1.T
        out[s, 1] = X2 @ X1.T
        out[s, 2] = X2 @ X2.T
    return out, n


def cv_spectrum(G):
    """G: (3, P, P) summed Grams. Symmetrised cvPCA spectrum, ranks 1..NR."""
    def one(A, C):
        w, U = np.linalg.eigh(A)
        U = U[:, ::-1][:, :NR]
        return np.einsum("pk,pq,qk->k", U, C, U) / (P - 1)
    return 0.5 * (one(G[0], G[1]) + one(G[2], G[1].T))


def spectrum_from(grams, members, splits):
    """grams: list per session of (S,3,P,P); members: session indices (with
    repeats); splits: which split indices to average."""
    lam = np.zeros(NR)
    for s in splits:
        G = np.zeros((3, P, P))
        for m in members:
            G += grams[m][s]
        lam += cv_spectrum(G)
    return lam / len(splits)


# --------------------------------------------------------------------------
# Forward model
# --------------------------------------------------------------------------

def single_pl(alpha, K):
    return np.arange(1, K + 1, dtype=np.float64) ** -alpha


def broken_pl(a1, a2, kb, K):
    k = np.arange(1, K + 1, dtype=np.float64)
    lam = k ** -a1
    hi = k > kb
    lam[hi] = kb ** (a2 - a1) * k[hi] ** -a2
    return lam


def model_library(specs, n_draws, seed):
    """specs: list of population spectra (1-D arrays). Returns mean and variance
    of the log cvPCA spectrum at 118 Gaussian images, ranks 1..NR."""
    rng = np.random.default_rng(seed)
    mu = np.zeros((len(specs), NR))
    m2 = np.zeros((len(specs), NR))
    for d in range(n_draws):
        Z = rng.standard_normal((P, EXACT_DIMS))
        Z -= Z.mean(0)
        C = np.eye(P) - 1.0 / P
        for j, lam in enumerate(specs):
            head = lam[:EXACT_DIMS]
            tail = lam[EXACT_DIMS:].sum()
            Zh = Z[:, :len(head)]
            G = (Zh * head) @ Zh.T + tail * C
            e = np.linalg.eigvalsh(G / (P - 1))[::-1][:NR]
            le = np.log(np.maximum(e, 1e-300))
            mu[j] += le
            m2[j] += le ** 2
    mu /= n_draws
    var = np.maximum(m2 / n_draws - mu ** 2, 0)
    return mu, var


def build_models(quick, seed):
    if quick:
        alphas = np.round(np.arange(0.6, 1.61, 0.05), 2)
        Ks = [1000, 10000]
        a1s, a2s, kbs = [0.5, 0.9], [1.0, 1.4], [10]
        n_draws = 6
    else:
        alphas = np.round(np.arange(0.50, 2.001, 0.02), 2)
        Ks = [300, 1000, 3000, 10000, 30000]
        a1s = np.round(np.arange(0.3, 1.301, 0.1), 1)
        a2s = np.round(np.arange(0.8, 2.001, 0.1), 1)
        kbs = [5, 10, 20]
        n_draws = 40
    sp = [(float(a), K) for K in Ks for a in alphas]
    bp = [(float(a1), float(a2), kb, K) for K in Ks for kb in kbs for a1 in a1s for a2 in a2s]
    key = hashlib.sha256(json.dumps([sp, bp, n_draws, seed, EXACT_DIMS]).encode()).hexdigest()[:16]
    f = CACHE / f"models_{key}.npz"
    if f.exists():
        z = np.load(f)
        return dict(single=sp, broken=bp, smu=z["smu"], svar=z["svar"], bmu=z["bmu"], bvar=z["bvar"], n_draws=n_draws)
    smu, svar = model_library([single_pl(a, K) for a, K in sp], n_draws, seed)
    bmu, bvar = model_library([broken_pl(a1, a2, kb, K) for a1, a2, kb, K in bp], n_draws, seed + 1)
    if not quick:
        f.parent.mkdir(parents=True, exist_ok=True)
        np.savez(f, smu=smu, svar=svar, bmu=bmu, bvar=bvar)
    return dict(single=sp, broken=bp, smu=smu, svar=svar, bmu=bmu, bvar=bvar, n_draws=n_draws)


def fit_range(lam, R):
    bad = np.flatnonzero(lam[:R] <= 0)
    return int(bad[0]) if len(bad) else R


def fit_models(lam, mu, w, R):
    """Weighted least squares with free additive constant over ranks 1..R.
    w: (models, NR) or (NR,) weights. Returns index of best model and errors."""
    r = fit_range(lam, R)
    y = np.log(lam[:r])
    m = mu[:, :r]
    ww = np.broadcast_to(w, mu.shape)[:, :r] if np.ndim(w) > 1 else np.broadcast_to(w[:r], m.shape)
    c = (ww * (y - m)).sum(1) / ww.sum(1)
    sse = (ww * (y - m - c[:, None]) ** 2).sum(1)
    return int(np.argmin(sse)), sse, r


def direct_slope(lam, lo, hi):
    k = np.arange(1, NR + 1)
    sel = (k >= lo) & (k <= hi) & (lam > 0)
    if sel.sum() < 3:
        return float("nan")
    return float(-np.polyfit(np.log(k[sel]), np.log(lam[sel]), 1)[0])


# --------------------------------------------------------------------------
# Per-area analysis
# --------------------------------------------------------------------------

def analyse_area(grams, mouse_of, models, rng, R=FIT_RANKS, n_boot=N_BOOT, n_cv=N_CV,
                 main_splits=None, do_cv=True):
    """grams: list per session (S,3,P,P); mouse_of: mouse id per session."""
    mice = sorted(set(mouse_of))
    by_mouse = {m: [i for i, x in enumerate(mouse_of) if x == m] for m in mice}
    S = grams[0].shape[0]
    main_splits = list(range(S)) if main_splits is None else main_splits
    all_sessions = list(range(len(grams)))
    lam = spectrum_from(grams, all_sessions, main_splits)

    boots = []
    for b in range(n_boot):
        pick = rng.choice(len(mice), len(mice), replace=True)
        members = [i for p in pick for i in by_mouse[mice[p]]]
        sp = rng.choice(S, min(2, S), replace=False)
        boots.append(spectrum_from(grams, members, sp))
    boots = np.array(boots)
    with np.errstate(invalid="ignore", divide="ignore"):
        lb = np.where(boots > 0, np.log(np.where(boots > 0, boots, 1)), np.nan)
    bvar = np.nanvar(lb, axis=0)
    bvar = np.where(np.isfinite(bvar), bvar, np.nanmax(bvar))

    smu, svar = models["smu"], models["svar"]
    w = 1.0 / (bvar[None, :] + svar)
    j, sse, r = fit_models(lam, smu, w, R)
    a_hat, K_hat = models["single"][j]
    a_boot, K_boot = [], []
    for bl in boots:
        jb, _, _ = fit_models(bl, smu, w, R)
        a_boot.append(models["single"][jb][0])
        K_boot.append(models["single"][jb][1])
    lo, hi = np.percentile(a_boot, [2.5, 97.5])

    wb = 1.0 / (bvar[None, :] + models["bvar"])
    jbk, bsse, _ = fit_models(lam, models["bmu"], wb, R)

    res = dict(
        n_mice=len(mice), n_sessions=len(grams),
        alpha_hat=a_hat, K_hat=K_hat, alpha_ci95=[float(lo), float(hi)],
        alpha_boot_median=float(np.median(a_boot)),
        K_boot_counts={str(k): int(np.sum(np.array(K_boot) == k)) for k in sorted(set(K_boot))},
        fit_ranks=r, single_wsse=float(sse[j]), broken_wsse=float(bsse[jbk]),
        broken_best=dict(zip(["alpha1", "alpha2", "k_break", "K"], models["broken"][jbk])),
        direct_slope_11_40=direct_slope(lam, 11, 40), direct_slope_1_10=direct_slope(lam, 1, 10),
        spectrum=[float(x) for x in lam],
    )
    if do_cv and len(mice) >= 4:
        wins, ratios = 0, []
        for c in range(n_cv):
            perm = rng.permutation(len(mice))
            hA, hB = perm[: len(mice) // 2], perm[len(mice) // 2:]
            mA = [i for p in hA for i in by_mouse[mice[p]]]
            mB = [i for p in hB for i in by_mouse[mice[p]]]
            sp = rng.choice(S, min(2, S), replace=False)
            lA, lB = spectrum_from(grams, mA, sp), spectrum_from(grams, mB, sp)
            one = np.ones(NR)
            js, _, rA = fit_models(lA, smu, one, R)
            jb2, _, _ = fit_models(lA, models["bmu"], one, R)
            rr = min(rA, fit_range(lB, R))
            y = np.log(lB[:rr])
            def err(mrow):
                d = y - mrow[:rr]
                return float(((d - d.mean()) ** 2).mean())
            es, eb = err(smu[js]), err(models["bmu"][jb2])
            wins += eb < es
            ratios.append(eb / es if es > 0 else np.nan)
        res["cv_broken_wins"] = int(wins)
        res["cv_median_error_ratio"] = float(np.nanmedian(ratios))
        res["broken_clearly_better"] = bool(wins >= 0.95 * n_cv and np.nanmedian(ratios) <= 0.9)
    else:
        res["broken_clearly_better"] = None
    res["verdict"] = verdict(res)
    return res


def verdict(res):
    lo, hi = res["alpha_ci95"]
    inside_T = T_BAND[0] <= lo and hi <= T_BAND[1]
    overlaps_B = lo <= B_BAND[1] and hi >= B_BAND[0]
    outside_T = hi < T_BAND[0] or lo > T_BAND[1]
    if outside_T:
        return "inconsistent"
    if inside_T and overlaps_B and not res.get("broken_clearly_better"):
        return "consistent"
    return "inconclusive"


# --------------------------------------------------------------------------
# Data: D1 two-photon
# --------------------------------------------------------------------------

def _get_json(url, path):
    import urllib.request
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        urllib.request.urlretrieve(url, path)
    return json.loads(path.read_text())


def head_info(url):
    import urllib.request
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return dict(etag=r.headers.get("ETag", "").strip('"'), size=int(r.headers.get("Content-Length", 0)),
                    last_modified=r.headers.get("Last-Modified"))


def d1_experiments():
    e = _get_json(S3 + "visual-coding-2p/ophys_experiments.json", CACHE / "d1" / "ophys_experiments.json")
    out = []
    for x in e:
        if x["stimulus_name"] != "three_session_B" or x["experiment_container"]["failed"]:
            continue
        area = x["targeted_structure"]["acronym"]
        if area not in AREAS:
            continue
        lines = [t["name"] for t in x["specimen"]["donor"]["transgenic_lines"]
                 if t["transgenic_line_type_name"] == "driver" and "Cre" in t["name"]]
        cre = lines[0] if lines else "unknown"
        out.append(dict(id=x["id"], area=area, mouse=x["specimen"]["donor_id"], depth=x["imaging_depth"],
                        cre=cre, excitatory=any(cre.startswith(g) for g in EXCITATORY)))
    return out


def _read_pandas_fixed(h, node):
    g = h[node]
    blocks = sorted(k for k in g.keys() if k.startswith("block") and k.endswith("_values"))
    cols, vals = [], []
    for b in blocks:
        items = [i.decode() if isinstance(i, bytes) else str(i) for i in g[b.replace("_values", "_items")][()]]
        cols += items
        vals.append(np.asarray(g[b][()]))
    return cols, np.concatenate(vals, axis=1)


def d1_load(exp_id):
    f = CACHE / "d1" / f"{exp_id}.npz"
    if f.exists():
        z = np.load(f, allow_pickle=False)
        return z["resp"], z["frame"], json.loads(str(z["info"]))
    import fsspec
    import h5py
    url = S3 + f"visual-coding-2p/ophys_experiment_analysis/{exp_id}_three_session_B_analysis.h5"
    info = head_info(url)
    fs = fsspec.filesystem("https")
    with fs.open(url, block_size=2 ** 20, cache_type="readahead") as fh, h5py.File(fh, "r") as h:
        cols, v = _read_pandas_fixed(h, "analysis/mean_sweep_response_ns")
        scol, sv = _read_pandas_fixed(h, "analysis/stim_table_ns")
    keep = [i for i, c in enumerate(cols) if c != "dx"]
    resp = v[:, keep].astype(np.float32)
    frame = sv[:, scol.index("frame")].astype(np.int64)
    f.parent.mkdir(parents=True, exist_ok=True)
    np.savez(f, resp=resp, frame=frame, info=json.dumps(dict(url=url, **info)))
    return resp, frame, dict(url=url, **info)


# --------------------------------------------------------------------------
# Data: D2 Neuropixels
# --------------------------------------------------------------------------

def d2_tables():
    import pandas as pd
    base = S3 + "visual-coding-neuropixels/ecephys-cache/"
    d = CACHE / "d2"
    d.mkdir(parents=True, exist_ok=True)
    import urllib.request
    for n in ["sessions.csv", "units.csv", "channels.csv", "probes.csv"]:
        if not (d / n).exists():
            urllib.request.urlretrieve(base + n, d / n)
    s = pd.read_csv(d / "sessions.csv")
    u = pd.read_csv(d / "units.csv")
    c = pd.read_csv(d / "channels.csv")
    p = pd.read_csv(d / "probes.csv")
    u = u.merge(c[["id", "ecephys_probe_id", "ecephys_structure_acronym"]].rename(columns={"id": "ecephys_channel_id"}),
                on="ecephys_channel_id").merge(
        p[["id", "ecephys_session_id"]].rename(columns={"id": "ecephys_probe_id"}), on="ecephys_probe_id")
    bo = s[s.session_type == "brain_observatory_1.1"]
    u = u[u.ecephys_session_id.isin(bo.id)]
    q = u[(u.quality == "good") & (u.isi_violations <= 0.5) & (u.amplitude_cutoff <= 0.1) & (u.presence_ratio >= 0.9)
          & u.ecephys_structure_acronym.isin(AREAS)]
    return bo, q


def d2_load(session_id, unit_ids):
    """Returns counts for windows [0,250) and [50,300) ms, frame, unit ids."""
    f = CACHE / "d2" / f"{session_id}.npz"
    if f.exists():
        z = np.load(f, allow_pickle=False)
        return z["c0"], z["c1"], z["frame"], z["uid"], json.loads(str(z["info"]))
    import fsspec
    import h5py
    url = S3 + f"visual-coding-neuropixels/ecephys-cache/session_{session_id}/session_{session_id}.nwb"
    info = head_info(url)
    fs = fsspec.filesystem("https")
    with fs.open(url, block_size=2 ** 22, cache_type="readahead") as fh, h5py.File(fh, "r") as h:
        ids = h["units/id"][()]
        idx = h["units/spike_times_index"][()]
        ns = h["intervals/natural_scenes_presentations"]
        start = ns["start_time"][()]
        frame = ns["frame"][()].astype(np.int64)
        want = set(int(x) for x in unit_ids)
        rows = [i for i, x in enumerate(ids) if int(x) in want]
        st = h["units/spike_times"]
        t0, t1 = start.min() - 1.0, start.max() + 1.0
        c0 = np.zeros((len(start), len(rows)), np.float32)
        c1 = np.zeros_like(c0)
        for j, i in enumerate(rows):
            a = int(idx[i - 1]) if i > 0 else 0
            b = int(idx[i])
            s = st[a:b]
            s = s[(s >= t0) & (s <= t1)]
            c0[:, j] = np.searchsorted(s, start + 0.25) - np.searchsorted(s, start)
            c1[:, j] = np.searchsorted(s, start + 0.30) - np.searchsorted(s, start + 0.05)
    uid = ids[rows].astype(np.int64)
    np.savez(f, c0=c0, c1=c1, frame=frame, uid=uid, info=json.dumps(dict(url=url, **info)))
    return c0, c1, frame, uid, dict(url=url, **info)


# --------------------------------------------------------------------------
# Session processing (worker)
# --------------------------------------------------------------------------

def complete(frame):
    imgs = set(np.unique(frame[frame >= 0]).tolist())
    counts = np.bincount(frame[frame >= 0], minlength=P)
    return imgs == set(range(P)) and counts.min() >= 2


def d1_worker(args):
    exp, seed = args
    resp, frame, info = d1_load(exp["id"])
    if not complete(frame):
        return exp, None, info
    rng = np.random.default_rng(seed)
    n = resp.shape[1]
    half = np.zeros(n, bool)
    half[rng.choice(n, n // 2, replace=False)] = True
    g_main, nm = session_grams(resp, frame, rng, N_SPLITS_MAIN)
    g_half, _ = session_grams(resp, frame, rng, N_SPLITS_VAR, neuron_mask=half)
    g_raw, _ = session_grams(resp, frame, rng, N_SPLITS_VAR, zscore=False)
    return exp, dict(main=g_main.astype(np.float32), half=g_half.astype(np.float32),
                     raw=g_raw.astype(np.float32), n=nm), info


def d2_worker(args):
    sid, units_by_area, seed = args
    allu = [u for us in units_by_area.values() for u in us]
    c0, c1, frame, uid, info = d2_load(sid, allu)
    if not complete(frame):
        return sid, None, info
    rng = np.random.default_rng(seed)
    out = {}
    for area, us in units_by_area.items():
        mask = np.isin(uid, us)
        if mask.sum() < 2:
            continue
        r0, r1 = c0[:, mask], c1[:, mask]
        n = r0.shape[1]
        half = np.zeros(n, bool)
        half[rng.choice(n, n // 2, replace=False)] = True
        g_main, nm = session_grams(r0, frame, rng, N_SPLITS_MAIN)
        g_half, _ = session_grams(r0, frame, rng, N_SPLITS_VAR, neuron_mask=half)
        g_raw, _ = session_grams(r0, frame, rng, N_SPLITS_VAR, zscore=False)
        g_win, _ = session_grams(r1, frame, rng, N_SPLITS_VAR)
        out[area] = dict(main=g_main.astype(np.float32), half=g_half.astype(np.float32),
                         raw=g_raw.astype(np.float32), win=g_win.astype(np.float32), n=nm)
    return sid, out, info


# --------------------------------------------------------------------------
# Synthetic data (quick mode and recovery test)
# --------------------------------------------------------------------------

def synthetic_sessions(alpha, K, n_mice, n_per, rng, reps=50, noise=1.5):
    """Each mouse: n_per neurons that are random mixtures of K stimulus features
    with variances k^-alpha, plus independent Gaussian trial noise."""
    lam = single_pl(alpha, K)
    Z = rng.standard_normal((P, K)) * np.sqrt(lam)
    frame = np.repeat(np.arange(P), reps)
    out = []
    for m in range(n_mice):
        W = rng.standard_normal((K, n_per)) / np.sqrt(K)
        S = Z @ W
        S /= S.std()
        resp = S[frame] + noise * rng.standard_normal((len(frame), n_per))
        out.append((resp.astype(np.float32), frame))
    return out


def run_synthetic(models, rng, alpha_true, K_true, n_mice, n_per, n_boot, n_cv):
    sess = synthetic_sessions(alpha_true, K_true, n_mice, n_per, rng)
    grams = [session_grams(r, f, rng, 4)[0] for r, f in sess]
    res = analyse_area(grams, list(range(n_mice)), models, rng, n_boot=n_boot, n_cv=n_cv)
    res.pop("spectrum")
    res["alpha_true"], res["K_true"] = alpha_true, K_true
    return res


# --------------------------------------------------------------------------

def versions():
    import importlib
    v = dict(python=platform.python_version(), platform=platform.platform())
    for m in ["numpy", "scipy", "h5py", "fsspec", "pandas", "aiohttp"]:
        try:
            v[m] = importlib.import_module(m).__version__
        except Exception:
            v[m] = None
    return v


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=20261125)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", choices=["d1", "d2"], default=None, help="debug: run one dataset")
    args = ap.parse_args()
    t_start = time.time()
    rng = np.random.default_rng(args.seed)
    print("building forward-model library ...", flush=True)
    models = build_models(args.quick, args.seed)
    print(f"  {len(models['single'])} single, {len(models['broken'])} broken models "
          f"({time.time() - t_start:.0f} s)", flush=True)

    if args.quick:
        res = run_synthetic(models, rng, 1.04, 3000, n_mice=8, n_per=150, n_boot=40, n_cv=10)
        ok = abs(res["alpha_hat"] - 1.04) <= 0.2
        print(json.dumps({k: res[k] for k in ["alpha_true", "alpha_hat", "K_hat", "alpha_ci95",
                                                "direct_slope_11_40", "verdict"]}, indent=1))
        print("quick synthetic recovery:", "PASS" if ok else "FAIL", f"({time.time() - t_start:.0f} s)")
        return 0 if ok else 1

    out = dict(seed=args.seed, versions=versions(), started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               bands=dict(claim=B_BAND, tolerance=T_BAND), n_boot=N_BOOT, n_cv=N_CV,
               model_draws=models["n_draws"])

    # Synthetic recovery (validation of the pipeline, not part of the verdict)
    print("synthetic recovery ...", flush=True)
    out["synthetic_recovery"] = [run_synthetic(models, rng, a, K, n_mice=20, n_per=200, n_boot=100, n_cv=20)
                                 for a, K in [(1.04, 3000), (1.04, 30000), (1.4, 3000), (0.8, 1000)]]
    for r in out["synthetic_recovery"]:
        print(f"  true a={r['alpha_true']} K={r['K_true']}: a_hat={r['alpha_hat']} CI={r['alpha_ci95']} "
              f"direct11-40={r['direct_slope_11_40']:.2f}", flush=True)

    results = {}
    files = []

    # ---- D1
    if args.only in (None, "d1"):
        exps = d1_experiments()
        print(f"D1: {len(exps)} session-B experiments in the six areas", flush=True)
        jobs = [(e, args.seed + int(e["id"]) % 100000) for e in exps]
        sess = []
        with ProcessPoolExecutor(args.workers) as ex:
            for i, (e, g, info) in enumerate(ex.map(d1_worker, jobs)):
                files.append(info)
                if g is None:
                    print("  excluded (incomplete images):", e["id"], flush=True)
                    continue
                sess.append((e, g))
                if (i + 1) % 25 == 0:
                    print(f"  {i + 1}/{len(jobs)} ({time.time() - t_start:.0f} s)", flush=True)
        d1 = {}
        for area in AREAS:
            for label, sel in [("main", lambda e: e["excitatory"]), ("all_cre", lambda e: True),
                               ("exploratory_L23_excitatory", lambda e: e["excitatory"] and e["depth"] < 250)]:
                if label != "main" and area != "VISp":
                    continue
                ss = [(e, g) for e, g in sess if e["area"] == area and sel(e)]
                if not ss:
                    continue
                grams = [g["main"] for e, g in ss]
                r = analyse_area(grams, [e["mouse"] for e, g in ss], models, rng, do_cv=(label == "main"))
                r["n_neurons"] = int(sum(g["n"] for e, g in ss))
                r["cre_lines"] = sorted(set(e["cre"] for e, g in ss))
                d1[f"{area}:{label}"] = r
                print(f"  D1 {area} {label}: N={r['n_neurons']} mice={r['n_mice']} a={r['alpha_hat']} "
                      f"CI={r['alpha_ci95']} K={r['K_hat']} verdict={r['verdict']}", flush=True)
            if area == "VISp":
                ss = [(e, g) for e, g in sess if e["area"] == area and e["excitatory"]]
                mo = [e["mouse"] for e, g in ss]
                for var in ["half", "raw"]:
                    r = analyse_area([g[var] for e, g in ss], mo, models, rng, do_cv=False)
                    d1[f"VISp:{var}"] = r
                    print(f"  D1 VISp {var}: a={r['alpha_hat']} CI={r['alpha_ci95']} verdict={r['verdict']}", flush=True)
                r = analyse_area([g["main"] for e, g in ss], mo, models, rng, R=60, do_cv=False)
                d1["VISp:ranks1-60"] = r
                print(f"  D1 VISp ranks1-60: a={r['alpha_hat']} CI={r['alpha_ci95']} verdict={r['verdict']}", flush=True)
        results["D1_two_photon"] = d1
        del sess

    # ---- D2
    if args.only in (None, "d2"):
        bo, q = d2_tables()
        jobs = []
        for sid in bo.id:
            uq = q[q.ecephys_session_id == sid]
            ub = {a: [int(x) for x in uq[uq.ecephys_structure_acronym == a].id] for a in AREAS}
            ub = {a: v for a, v in ub.items() if v}
            if ub:
                jobs.append((int(sid), ub, args.seed + int(sid) % 100000))
        print(f"D2: {len(jobs)} Neuropixels sessions", flush=True)
        per_area = {a: [] for a in AREAS}
        with ProcessPoolExecutor(min(args.workers, 6)) as ex:
            for i, (sid, g, info) in enumerate(ex.map(d2_worker, jobs)):
                files.append(info)
                if g is None:
                    print("  excluded (incomplete images):", sid, flush=True)
                    continue
                for a, ga in g.items():
                    per_area[a].append((sid, ga))
                print(f"  {i + 1}/{len(jobs)} ({time.time() - t_start:.0f} s)", flush=True)
        d2 = {}
        for area in AREAS:
            ss = per_area[area]
            if not ss:
                continue
            mo = [s for s, g in ss]
            r = analyse_area([g["main"] for s, g in ss], mo, models, rng)
            r["n_neurons"] = int(sum(g["n"] for s, g in ss))
            d2[f"{area}:main"] = r
            print(f"  D2 {area}: N={r['n_neurons']} mice={r['n_mice']} a={r['alpha_hat']} CI={r['alpha_ci95']} "
                  f"K={r['K_hat']} verdict={r['verdict']}", flush=True)
            if area == "VISp":
                for var in ["half", "raw", "win"]:
                    r = analyse_area([g[var] for s, g in ss], mo, models, rng, do_cv=False)
                    d2[f"VISp:{var}"] = r
                    print(f"  D2 VISp {var}: a={r['alpha_hat']} CI={r['alpha_ci95']} verdict={r['verdict']}", flush=True)
                r = analyse_area([g["main"] for s, g in ss], mo, models, rng, R=60, do_cv=False)
                d2["VISp:ranks1-60"] = r
                print(f"  D2 VISp ranks1-60: a={r['alpha_hat']} CI={r['alpha_ci95']} verdict={r['verdict']}", flush=True)
        results["D2_neuropixels"] = d2

    out["results"] = results
    out["files"] = files

    if args.only is None:
        d1, d2 = results["D1_two_photon"], results["D2_neuropixels"]
        v1, v2 = d1["VISp:main"]["verdict"], d2["VISp:main"]["verdict"]
        flags = {}
        leaves = False   # a planned robustness interval is not inside T
        flips = False    # a planned robustness check flips consistent <-> inconsistent
        planned = {"VISp:half", "VISp:raw", "VISp:ranks1-60", "VISp:all_cre", "VISp:win"}
        for name, d, base in [("D1", d1, v1), ("D2", d2, v2)]:
            for k, r in d.items():
                if k not in planned:
                    continue
                lo, hi = r["alpha_ci95"]
                not_inside = lo < T_BAND[0] or hi > T_BAND[1]
                fl = {base, r["verdict"]} == {"consistent", "inconsistent"}
                flags[f"{name} {k}"] = dict(verdict=r["verdict"], ci=r["alpha_ci95"], alpha_hat=r["alpha_hat"],
                                            interval_not_inside_T=bool(not_inside), flips=bool(fl))
                if k == "VISp:half" and abs(r["alpha_hat"] - d["VISp:main"]["alpha_hat"]) > 0.10:
                    flags[f"{name} {k}"]["large_population_assumption_not_met"] = True
                leaves |= not_inside
                flips |= fl
        both_broken = bool(d1["VISp:main"].get("broken_clearly_better")) and bool(d2["VISp:main"].get("broken_clearly_better"))
        if both_broken:
            status = "Contested"
        elif v1 == "consistent" and v2 == "consistent" and not leaves:
            status = "Partial (consistent)"
        elif v1 == "inconsistent" and v2 == "inconsistent" and not flips:
            status = "Refuted"
        elif {v1, v2} == {"consistent", "inconsistent"} or flips:
            status = "Contested"
        else:
            status = "Partial (inconclusive)"
        out["decision"] = dict(D1_V1=v1, D2_V1=v2, robustness=flags, robustness_interval_not_inside_T=leaves, robustness_flip=flips,
                               broken_clearly_better_in_both=both_broken, status=status)
        print("DECISION:", json.dumps(out["decision"], indent=1))
    out["runtime_s"] = round(time.time() - t_start, 1)
    name = "check-output.json" if args.only is None else f"check-output-{args.only}.json"
    (HERE / name).write_text(json.dumps(out, indent=1, default=float) + "\n")
    if args.only is None:
        return 0 if out["decision"]["status"] == "Partial (consistent)" else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
