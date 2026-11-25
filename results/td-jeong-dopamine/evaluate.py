"""Score a TD model on Tests 1 to 11 with the frozen pass rules (METHODS.md)."""
from __future__ import annotations

import numpy as np

import jeong as J
import tests as T

NSEEDS = 20


# ------------------------------------------------------------------ measures on any source
def measures_for(get, sess, groups):
    M = {}
    if "T12" in groups:
        M["T12"] = {s: T.measure_randomrewards(get, sess[s], s) for s in J.RR8}
    for iC, key in ((0, "T3"), (1, "T4"), (2, "T6")):
        if key in groups:
            M[key] = {s: T.measure_cumsum(get, sess[s], s, iC) for s in J.MAIN7}
    for key, f in (("T5", T.measure_T5), ("T7", T.measure_T7), ("T8", T.measure_T8), ("T9", T.measure_T9)):
        if key in groups:
            M[key] = {s: f(get, sess[s], s) for s in J.MAIN7}
    if "SEQ" in groups:
        M["SEQ_CTRL"] = {s: T.measure_seq(get, sess[s], s) for s in J.SEQ_CTRL}
        M["SEQ_INH"] = {s: T.measure_seq(get, sess[s], s) for s in J.SEQ_INH}
    return M


def scale_ref(get, sess, subj):
    """Reward response used to put model and data in the same units."""
    vals = []
    if subj in J.MAIN7:
        pav = [(d, n) for d, n, p in sess if p == "Pavlovian" and d < J.LEARNEDDAY[J.MAIN7.index(subj)]][-3:]
        for d, n in pav:
            s = get(subj, n)
            e = J.eventfrompavlovian(s)
            nn = min(len(e["CStime"]), len(e["firstfxlicktime"]))
            first, cst = e["firstfxlicktime"][:nn], e["CStime"][:nn]
            tr, rw = J.alignsignal2event(s["t"], s["sig"], first, [-2000, 5000], 10)
            tc, cs = J.alignsignal2event(s["t"], s["sig"], cst, [-2000, 5000], 10)
            v = J.aucsignal(rw, tr, [0, 500]) - J.aucsignal(cs, tc, [-500, 0])
            vals.append(v[~np.isnan(first)])
    elif "stGtACR" in subj:
        pav = [(d, n) for d, n, p in sess if p == "Pavlovian"][:1]
        for d, n in pav:
            s = get(subj, n)
            et, ei = s["eventtime"], s["eventindex"]
            fl = T.first_after(et, ei, J.LICK, J.FXRW, 5000)
            cs1 = et[ei == J.CS1][0::2]
            b = np.nanmean(T.calc_auc_s(s["sig"], s["t"], cs1 - 500, [0, 500]))
            vals.append(T.calc_auc_s(s["sig"], s["t"], fl, [0, 500]) - b)
    else:  # M8: RandomRewards only
        rr = [(d, n) for d, n, p in sess if p == "RandomRewards"][-3:]
        for d, n in rr:
            s = get(subj, n)
            et, ei = s["eventtime"], s["eventindex"]
            fl = T.first_after(et, ei, J.LICK, J.BGD, 5000)
            vals.append(J.aucsignal2event(s["t"], s["sig"], fl, [0, 500])
                        - J.aucsignal2event(s["t"], s["sig"], fl, [-2000, -1500]))
    v = np.concatenate(vals) if vals else np.array([np.nan])
    return float(np.nanmean(v))


# ------------------------------------------------------------------ noise and scaling
def _walk(obj, fn, path=()):
    """Apply fn(path, array) to every 'da_' array in nested dict/list measures."""
    if isinstance(obj, dict):
        return {k: (fn(path + (k,), v) if isinstance(k, str) and k.startswith("da_") and isinstance(v, (np.ndarray, list, float, np.floating))
                    else _walk(v, fn, path + (k,))) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_walk(v, fn, path + (i,)) for i, v in enumerate(obj)]
    return obj


def _resid_sd(x):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if len(x) < 5:
        return float(np.std(x)) if len(x) > 1 else 0.0
    from scipy.ndimage import median_filter
    return float(np.std(x - median_filter(x, size=21, mode="nearest")))


def noise_sds(Mdata):
    """SD per (group, animal, measure key) from the real data; sessions pooled."""
    sds = {}
    for g, byanimal in Mdata.items():
        for subj, m in byanimal.items():
            pool = {}

            def grab(path, v):
                key = path[-1]
                if isinstance(v, list):
                    for a in v:
                        pool.setdefault(key, []).append(np.atleast_1d(np.asarray(a, float)))
                else:
                    pool.setdefault(key, []).append(np.atleast_1d(np.asarray(v, float)))
                return v
            _walk(m, grab)
            for key, arrs in pool.items():
                sds[(g, subj, key)] = _resid_sd(np.concatenate(arrs))
    return sds


def add_noise(Mmodel, scales, sds, seed):
    rng = np.random.default_rng(seed)
    out = {}
    for g, byanimal in Mmodel.items():
        out[g] = {}
        for subj, m in byanimal.items():
            sc = scales[subj]

            def f(path, v, g=g, subj=subj, sc=sc):
                sd = sds.get((g, subj, path[-1]), 0.0)
                if isinstance(v, list):
                    return [np.asarray(a, float) * sc + rng.normal(0, sd, np.shape(a)) for a in v]
                a = np.asarray(v, float)
                return a * sc + rng.normal(0, sd, a.shape)
            out[g][subj] = _walk(m, f)
    return out


# ------------------------------------------------------------------ stats and pass rules
def stats_of(M):
    st = {}
    if "T12" in M:
        st["T1"], st["T2"] = T.stat_T1(M["T12"]), T.stat_T2(M["T12"])
    for key, f in (("T3", T.stat_T3), ("T4", T.stat_T4), ("T5", T.stat_T5), ("T6", T.stat_T6),
                   ("T7", T.stat_T7), ("T8", T.stat_T8), ("T9", T.stat_T9)):
        if key in M:
            st[key] = f(M[key])
    if "SEQ_CTRL" in M:
        st["T10"], st["T11"] = T.stat_T10(M["SEQ_CTRL"]), T.stat_T11(M["SEQ_INH"])
    return st


def sig_same(r, sign):
    return bool(np.sign(r["t"]) == sign and r["p"] < 0.05)


def passes(st, data_st):
    """Frozen pass rules. data_st gives the data CIs for Tests 9 and 10."""
    P = {}
    if "T1" in st:
        P["T1"] = sig_same(st["T1"], +1)
        P["T2"] = sig_same(st["T2"], +1)
    if "T3" in st:
        P["T3"] = sig_same(st["T3"]["abrupt"], +1) and sig_same(st["T3"]["change"], +1)  # deviation 2
    if "T4" in st:
        r = st["T4"]["daic_da"]
        P["T4"] = bool(r["p"] >= 0.05 or r["t"] <= 0)
    if "T5" in st:
        sgn = np.sign(data_st["T5"]["paired"]["t"])
        od = st["T5"]["old_delay"]
        P["T5"] = sig_same(st["T5"]["paired"], sgn) and bool(od["t"] >= 0 or od["p"] >= 0.05)
    if "T6" in st:
        P["T6"] = sig_same(st["T6"]["abrupt"], +1) and sig_same(st["T6"]["mean_da"], +1)
    if "T7" in st:
        P["T7"] = sig_same(st["T7"], -1)
    if "T8" in st:
        P["T8"] = sig_same(st["T8"]["cue_ratio"], +1) and sig_same(st["T8"]["reward_ratio"], +1)
    if "T9" in st:
        lo, hi = data_st["T9"]["ci"]
        P["T9"] = bool(lo <= st["T9"]["mean"] <= hi)
    if "T10" in st:
        lo, hi = data_st["T10"]["ci"]
        P["T10"] = bool(lo <= st["T10"]["mean"] <= hi)
        P["T11"] = bool(st["T11"]["all_positive"])
    return P


def median_stats(stats_list):
    """Median over noise seeds of every t, p and mean (recursively)."""
    def med(xs):
        x0 = xs[0]
        if isinstance(x0, dict):
            return {k: med([x[k] for x in xs]) for k in x0}
        if isinstance(x0, (int, float, np.floating, np.integer)) and not isinstance(x0, bool):
            return float(np.nanmedian([float(x) for x in xs]))
        if isinstance(x0, bool):
            return bool(np.mean(xs) >= 0.5)
        if isinstance(x0, list) and x0 and isinstance(x0[0], (float, int, np.floating)):
            return [float(v) for v in np.nanmedian(np.array(xs, float), axis=0)]
        return x0
    return med(stats_list)


def margins(st, data_st):
    """Signed distance of each scored criterion from its pass boundary (>0 passes),
    in t units (or CI half-widths for Tests 9, 10). Used only to rank failing models."""
    from scipy import stats as SS

    def tc(n):
        return float(SS.t.ppf(0.975, max(n, 2) - 1))
    m = {}

    def sig(name, r, sign):
        m[name] = sign * r["t"] - tc(r["n"]) if np.isfinite(r["t"]) else -np.inf
    if "T1" in st:
        sig("T1", st["T1"], 1); sig("T2", st["T2"], 1)
    if "T3" in st:
        sig("T3a", st["T3"]["abrupt"], 1); sig("T3c", st["T3"]["change"], 1)
    if "T4" in st:
        r = st["T4"]["daic_da"]
        m["T4"] = tc(r["n"]) - r["t"] if np.isfinite(r["t"]) else -np.inf
    if "T5" in st:
        sig("T5p", st["T5"]["paired"], np.sign(data_st["T5"]["paired"]["t"]))
        r = st["T5"]["old_delay"]
        m["T5o"] = r["t"] + tc(r["n"]) if np.isfinite(r["t"]) else -np.inf
    if "T6" in st:
        sig("T6a", st["T6"]["abrupt"], 1); sig("T6m", st["T6"]["mean_da"], 1)
    if "T7" in st:
        sig("T7", st["T7"], -1)
    if "T8" in st:
        sig("T8c", st["T8"]["cue_ratio"], 1); sig("T8r", st["T8"]["reward_ratio"], 1)
    for k in ("T9", "T10"):
        if k in st:
            lo, hi = data_st[k]["ci"]
            h = (hi - lo) / 2
            x = st[k]["mean"]
            m[k] = (min(x - lo, hi - x) / h) if np.isfinite(x) else -np.inf
    if "T11" in st:
        v = np.array(st["T11"]["last_cs1_over_rw"], float)
        m["T11"] = float(np.nanmin(v) / (np.nanstd(v) + 1e-9)) if np.isfinite(v).any() else -np.inf
    return m
