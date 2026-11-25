"""Jeong et al. 2022 measures, ported from namboodirilab/ANCCR analysis code.

Everything runs in milliseconds, as the MATLAB code does, so window and bin
rounding match. A session is a dict with eventindex, eventtime (ms),
nonsolenoidflag, t (dF/F timestamps, ms) and sig (dF/F). The same functions
are applied to real dF/F and to a model's simulated dF/F.
"""
from __future__ import annotations

import os
import re

import numpy as np

# ----------------------------------------------------------------- animals
MAIN7 = ["HJ-FP-M2", "HJ-FP-M3", "HJ-FP-M4", "HJ-FP-F1", "HJ-FP-F2", "HJ-FP-M6", "HJ-FP-M7"]
RR8 = MAIN7 + ["HJ-FP-M8"]
SEQ_CTRL = ["HJ-FP-datWT-stGtACR-F1"] + ["HJ-FP-WT-stGtACR-" + i for i in ["F1", "F2", "F3", "M1", "M2", "M3"]]
SEQ_INH = ["HJ-FP-datHT-stGtACR-" + i for i in ["M2", "M3", "M4", "M5", "M6", "M7"]]

# fig4/cumsum_pavlovian.m, rows in MAIN7 order, columns = condition 1..4
STARTDAY = np.array([[16, 29, 70, 38], [12, 24, 61, 28], [13, 32, 66, 38], [7, 19, 42, 23],
                     [13, 24, 46, 28], [7, 19, 47, 23], [6, 20, 43, 24]])
ENDDAY = np.array([[28, 32, 72, 43], [23, 27, 63, 38], [31, 37, 67, 47], [18, 22, 43, 27],
                   [23, 27, 47, 33], [18, 22, 50, 30], [19, 23, 44, 29]])
LEARNEDDAY = np.array([22, 17, 25, 12, 17, 12, 14])
RATIO = [1.5, 1.9, 2.5, 2]
NDAYBEFORE = [0, 3, 1, 1]
# fig4/cumsum_bgdrw_vs_extinction.m
BGD_START = np.array([[33, 70], [28, 61], [38, 66], [23, 42], [28, 46], [23, 47], [24, 43]])
BGD_END = np.array([[43, 72], [38, 63], [47, 67], [27, 43], [33, 47], [30, 50], [29, 44]])
# figS10/omission.m
OMIT_DAY = np.array([[29, 24, 32, 19, 24, 19, 20], [70, 61, 66, 42, 46, 47, 43]])
# fig5/intermediate_rsp_trialless.m
POIS_LEARNED = np.array([77, 67, 70, 46, 51, 53, 48])
POIS_END = np.array([78, 68, 72, 49, 54, 55, 50])
# fig6/backpropagation_pavlovian.m
BP_START = np.array([16, 12, 13, 7, 13, 7, 6])
BP_END = np.array([28, 23, 31, 18, 23, 18, 19])

LICK, BGD, FXRW, TRIALEND, CS1, CS2, SESEND = 5, 7, 10, 14, 15, 16, 0


# ----------------------------------------------------------------- loading
def day_of(name: str) -> int:
    return int(re.match(r"Day(\d+)", name).group(1))


def list_sessions(cache: str, subj: str):
    """[(day, name, protocol)] sorted by day; protocol from session_description."""
    out = []
    for fn in os.listdir(os.path.join(cache, subj)):
        if fn.endswith(".npz"):
            name = fn[:-4]
            with np.load(os.path.join(cache, subj, fn)) as z:
                prot = str(z["session_description"])
            out.append((day_of(name), name, prot))
    return sorted(out)


def load(cache: str, subj: str, name: str) -> dict:
    with np.load(os.path.join(cache, subj, name + ".npz")) as z:
        d = {k: z[k] for k in z.files}
    return dict(eventindex=d["eventindex"].astype(int), eventtime=d["eventtime"] * 1000.0,
                nonsolenoidflag=d["nonsolenoidflag"].astype(int), t=d["dff_t"] * 1000.0,
                sig=d["dff"].astype(np.float64), protocol=str(d["session_description"]),
                desc=d["desc"], code=d["code"])


# ----------------------------------------------------------------- MATLAB ports
def gauss_kernel(n: int, sigma: float) -> np.ndarray:
    x = np.arange(n) - (n - 1) / 2
    h = np.exp(-x ** 2 / (2 * sigma ** 2))
    return h / h.sum()


def conv_same(a: np.ndarray, k: np.ndarray) -> np.ndarray:
    full = np.convolve(a, k)
    s = len(k) // 2
    return full[s:s + len(a)]


def conv_valid(a: np.ndarray, k: np.ndarray) -> np.ndarray:
    return np.convolve(a, k, mode="valid")


def nearest_idx(t: np.ndarray, x: np.ndarray) -> np.ndarray:
    """MATLAB min(abs(t - x)) index (first index on ties); NaNs in t never win; -1 for NaN x."""
    x = np.asarray(x, float)
    good = np.nonzero(~np.isnan(t))[0]
    tg = t[good]
    if np.all(np.diff(tg) >= 0):
        i = np.searchsorted(tg, np.where(np.isnan(x), 0, x))
        lo = np.clip(i - 1, 0, len(tg) - 1)
        hi = np.clip(i, 0, len(tg) - 1)
        pick = np.where(np.abs(tg[lo] - x) <= np.abs(tg[hi] - x), lo, hi)
        out = good[pick]
        out[np.isnan(x)] = -1
        return out
    out = np.zeros(len(x), dtype=int)
    for j, v in enumerate(x):
        out[j] = -1 if np.isnan(v) else int(np.nanargmin(np.abs(t - v)))
    return out


def aucsignal2event(t, sig, ev, window, out=None):
    ev = np.asarray(ev, dtype=float)
    out = np.zeros(len(ev), bool) if out is None else np.asarray(out, bool).copy()
    dt = np.nanmean(np.diff(t))
    w0, w1 = int(np.round(window[0] / dt)), int(np.round(window[1] / dt))
    out |= np.isnan(ev)
    idx = nearest_idx(t, ev)
    res = np.full(len(ev), np.nan)
    for j in range(len(ev)):
        if out[j]:
            continue
        i = idx[j]
        if i + w1 > len(t) - 1 or i + w0 < 0:
            continue
        res[j] = np.sum(sig[i + w0:i + w1 + 1])
    return res


def alignsignal2event(t, sig, ev, window, resolution):
    ev = np.asarray(ev, dtype=float)
    dt = np.round(np.nanmean(np.diff(t)) * 10) / 10
    w0, w1 = int(np.round(window[0] / dt)), int(np.round(window[1] / dt))
    time = np.arange(w0, w1 + 1) * dt
    idx = nearest_idx(t, ev)
    n = len(sig)
    cols = idx[:, None] + np.arange(w0, w1 + 1)[None, :]
    valid = (cols >= 0) & (cols < n) & (idx[:, None] >= 0)
    al = np.where(valid, sig[np.clip(cols, 0, n - 1)], np.nan)
    if resolution > 0 and len(al):
        k = gauss_kernel(5 * resolution, resolution)
        from scipy.signal import fftconvolve  # noqa: F401
        full = np.apply_along_axis(lambda r: np.convolve(r, k), 1, al)
        st = len(k) // 2
        al = full[:, st:st + al.shape[1]]
    return time, al


def aucsignal(signal, time, window, out=None):
    binsize = np.nanmean(np.diff(time))
    m = (time >= window[0]) & (time <= window[1])
    res = np.nansum(signal[:, m], axis=1) * binsize
    if out is not None:
        res = res.astype(float)
        res[np.asarray(out, bool)] = np.nan
    return res


def firsttimeafterevent(target, ref, window):
    target = np.sort(np.asarray(target, float))
    ref = np.asarray(ref, float)
    res = np.full(len(ref), np.nan)
    if window is None:
        for j, x in enumerate(ref):
            i = np.searchsorted(target, x, side="left")
            if i < len(target):
                res[j] = target[i]
        return res
    w = np.asarray(window, float)
    if w.ndim == 0 or len(w) == 1:
        w = ref + float(np.ravel(w)[0])
    for j, (x, y) in enumerate(zip(ref, w)):
        i = np.searchsorted(target, x, side="right")
        if i < len(target) and target[i] < y:
            res[j] = target[i]
    return res


def numevent(ev, ref, window, normalize):
    ev = np.asarray(ev, float)
    n = np.array([np.sum((ev >= x + window[0]) & (ev <= x + window[1])) for x in ref], float)
    if normalize == 1:
        n = n / ((window[1] - window[0]) / 1000)
    return n


def point2line(pt, v1, v2):
    a = np.array([v1[0] - v2[0], v1[1] - v2[1], 0.0])
    b = np.column_stack([pt[:, 0] - v2[0], pt[:, 1] - v2[1], np.zeros(len(pt))])
    return np.linalg.norm(np.cross(a, b), axis=1) / np.linalg.norm(a)


def eventfrompavlovian(s, trialless=False):
    et, ei, ns = s["eventtime"], s["eventindex"], s["nonsolenoidflag"]
    sessionend = et[ei == SESEND]
    sessionend = sessionend if len(sessionend) else np.array([et[-1] + 5000])
    licktime = et[ei == LICK]
    csmask = np.isin(ei, [CS1, CS2])
    CS = np.where(ei[csmask] == CS1, 1, 2)
    CStime = et[csmask]
    fxrwtime = et[ei == FXRW]
    fxreward = 1 - ns[ei == FXRW]
    bgdrwtime = et[ei == BGD]
    keep = CStime <= fxrwtime[-1]
    CS, CStime = CS[keep], CStime[keep]
    n = min(len(CStime), len(fxrwtime))
    CSduration = np.nanmean(fxrwtime[:n] - CStime[:n]) - 1000
    trialend = et[ei == TRIALEND]
    if not trialless:
        if len(trialend) and np.nanmean(trialend[:len(fxrwtime)] - fxrwtime[:len(trialend)]) >= 2000:
            first = firsttimeafterevent(licktime, fxrwtime, trialend[:len(fxrwtime)])
        else:
            first = firsttimeafterevent(licktime, fxrwtime, np.r_[CStime[1:] - 1000, sessionend[:1]])
    else:
        first = firsttimeafterevent(licktime, fxrwtime, None)
    return dict(sessionend=sessionend[0], licktime=licktime, CS=CS, CStime=CStime, CSduration=CSduration,
                fxrwtime=fxrwtime, fxreward=fxreward, bgdrwtime=bgdrwtime, firstfxlicktime=first)
