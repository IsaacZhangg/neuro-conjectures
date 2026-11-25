"""Jeong et al. 2022 Tests 1 to 11: per-trial measures and group statistics.

measure_* functions read sessions through `get(subj, name)` (data or model, same
format) and return per-trial arrays. Dopamine arrays are under keys starting
with "da_"; anything else (licks, intervals) comes from the event log only.
stat_* functions take {subj: measures} and return the group statistic.
Ported from namboodirilab/ANCCR analysis/{fig3,fig4,fig5,fig6,figS10}.
"""
from __future__ import annotations

import numpy as np
from scipy import optimize, stats

import jeong as J

# ------------------------------------------------------------------ helpers


def sessions_in(sess, lo, hi, protocol=None):
    return [(d, n) for d, n, p in sess if lo <= d <= hi and (protocol is None or p == protocol)]


def ttest1(x, mu=0.0):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    r = stats.ttest_1samp(x, mu)
    return dict(t=float(r.statistic), p=float(r.pvalue), n=len(x), mean=float(np.mean(x)),
                ci=[float(v) for v in stats.t.interval(0.95, len(x) - 1, loc=np.mean(x), scale=stats.sem(x))])


def ttest_paired(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    r = stats.ttest_rel(a, b)
    d = a - b
    return dict(t=float(r.statistic), p=float(r.pvalue), n=len(a), mean=float(np.mean(d)),
                ci=[float(v) for v in stats.t.interval(0.95, len(d) - 1, loc=np.mean(d), scale=stats.sem(d))])


def matlab_corr(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 3:
        return np.nan
    return float(np.corrcoef(x[m], y[m])[0, 1])


# ------------------------------------------------------------------ Tests 1, 2


def measure_randomrewards(get, sess, subj):
    """fig3/corr_randomrewards.m. 12 s IRI sessions only (stop at the first 30 s one)."""
    out = []
    for d, n, p in sess:
        if p != "RandomRewards":
            continue
        s = get(subj, n)
        et, ei = s["eventtime"], s["eventindex"]
        bgd = et[ei == J.BGD]
        if len(bgd) < 5:
            continue
        if np.mean(np.diff(bgd)) > 20000:  # first 30 s IRI session ends the set
            break
        lick = et[ei == J.LICK]
        send = et[ei == J.SESEND]
        send = send[:1] if len(send) else np.array([et[-1] + 5000])
        iri = np.r_[np.nan, np.diff(bgd)]
        first = J.firsttimeafterevent(lick, bgd, np.r_[bgd[1:], send])
        excl = (iri < 3000) | np.isnan(first)
        excl[0] = False
        a = J.aucsignal2event(s["t"], s["sig"], first, [-500, 1000], excl)
        b = J.aucsignal2event(s["t"], s["sig"], first, [-2000, -500], excl)
        out.append(dict(day=d, da_auc=(a - b)[~excl], iri=iri[~excl], idx=np.nonzero(~excl)[0], nrw=len(bgd)))
    return out


def stat_T1(M):
    r = []
    for subj, ses in M.items():
        n0, num, auc = 0, [], []
        for s in ses:
            num.append(s["idx"] + 1 + n0)
            auc.append(s["da_auc"])
            n0 += s["nrw"]
        r.append(matlab_corr(np.concatenate(num), np.concatenate(auc)))
    return ttest1(r) | dict(per_animal=r)


def stat_T2(M):
    r = [np.nanmean([matlab_corr(s["iri"], s["da_auc"]) for s in ses]) for ses in M.values()]
    return ttest1(r) | dict(per_animal=[float(x) for x in r])


# ------------------------------------------------------------------ Tests 3, 4, 6 (cumsum_pavlovian.m)


def cs_auc(s, ev, base=(-1000, 0), win=(0, 2000), div=2.0, res=10, span=(-2000, 10000)):
    time, al = J.alignsignal2event(s["t"], s["sig"], ev, list(span), res)
    return J.aucsignal(al, time, list(win)) / div - J.aucsignal(al, time, list(base))


def measure_cumsum(get, sess, subj, iC):
    """iC = 0 initial acquisition, 1 delay extension to 8 s cue, 2 extinction (Jeong iC 1..3)."""
    iM = J.MAIN7.index(subj)
    pav = [(d, n) for d, n, p in sess if p == "Pavlovian"]
    lo = J.STARTDAY[iM, iC] - J.NDAYBEFORE[iC]
    hi = J.ENDDAY[iM, iC]
    first_pav = get(subj, pav[0][1])
    e0 = J.eventfrompavlovian(first_pav)
    csrw = np.unique(e0["CS"][: len(e0["fxreward"])][e0["fxreward"][: len(e0["CS"])] == 1])[0]
    days, lick, base, auc = [], [], [], []
    for d, n in pav:
        if not (lo <= d <= hi):
            continue
        s = get(subj, n)
        e = J.eventfrompavlovian(s)
        cst = e["CStime"][e["CS"] == csrw]
        a = cs_auc(s, cst)
        if iC == 1:
            lt = e["licktime"]
            lc = np.array([np.nanmedian(lt[(lt >= x) & (lt <= x + e["CSduration"] + 1000)] - x)
                           if np.any((lt >= x) & (lt <= x + e["CSduration"] + 1000)) else np.nan for x in cst])
            keep = ~np.isnan(lc)
            lick.append(lc[keep]); a = a[keep]; base.append(np.zeros(keep.sum()))
        else:
            lick.append(J.numevent(e["licktime"], cst, [0, e["CSduration"] + 1000], 1))
            base.append(J.numevent(e["licktime"], cst, [-1000, 0], 1))
        auc.append(a)
        days.append(np.full(len(a), d))
    return dict(days=np.concatenate(days), lick=np.concatenate(lick), lickbase=np.concatenate(base),
                da_auc=np.concatenate(auc), iC=iC, iM=iM)


def weibull(x, A, L, S, b):
    return A * (1 - 2.0 ** (-(x / L) ** S)) + b


def fit_aic(y, lower, upper, nbefore, seed=7):
    n = len(y)
    x = np.arange(1, n + 1, dtype=float)
    rng = np.random.default_rng(seed)
    best = np.inf
    lo, hi = np.array(lower, float), np.array(upper, float)
    for k in range(20):
        p0 = lo + rng.random(4) * (hi - lo)
        if k == 0:
            p0 = np.clip([np.nanmax(y) - np.nanmin(y), nbefore + 1, 2.0, np.nanmin(y)], lo + 1e-9, hi - 1e-9)
        try:
            p, _ = optimize.curve_fit(weibull, x, y, p0=p0, bounds=(lo, hi), maxfev=20000)
            rss = np.sum((y - weibull(x, *p)) ** 2)
            best = min(best, rss)
        except Exception:
            continue
    aic_w = 2 * 4 + n * np.log(best)
    aic_c = 2 * 1 + n * np.log(np.sum((y - y.mean()) ** 2))
    return aic_c - aic_w  # diff(daic,[],2): constant minus Weibull


def cumsum_stats(m):
    iC, iM = m["iC"], m["iM"]
    days = m["days"]
    inbefore = days < (J.LEARNEDDAY[iM] if iC == 0 else J.STARTDAY[iM, iC])
    nbefore = int(inbefore.sum())
    ntot = int(round(nbefore * J.RATIO[iC]))
    ntot = min(ntot, len(days))
    lickcum_cs = np.cumsum(m["lick"])
    lickcum = lickcum_cs if iC == 1 else lickcum_cs - np.cumsum(m["lickbase"])
    auccum = np.cumsum(m["da_auc"])
    k = J.gauss_kernel(10, 2)
    out = {}
    for key, cum in (("lick", lickcum), ("da", auccum)):
        c = J.conv_valid(cum[:ntot], k)
        c = c / c[-1]
        xs = np.arange(5, ntot - 5 + 1)
        c = c[: len(xs)] if len(c) >= len(xs) else np.r_[c, np.full(len(xs) - len(c), c[-1])]
        d = J.point2line(np.column_stack([xs, c]), [5, c[0]], [ntot - 5, 1])
        out["abrupt_" + key] = float(np.max(d))
        out["change_" + key] = int(np.argmax(d)) + 1
    if iC == 1:
        out["daic_lick"] = fit_aic(m["lick"][:ntot], [0, nbefore - 10, 0, 0], [9000, ntot, 100, 4000], nbefore)
        out["daic_da"] = fit_aic(m["da_auc"][:ntot], [-2, nbefore - 10, 0, -10], [50, ntot, 10, 10], nbefore)
    if iC == 2:
        out["mean_da_ext"] = float(np.mean(m["da_auc"][nbefore:ntot]))
    return out


def stat_T3(M):
    r = [cumsum_stats(m) for m in M.values()]
    return dict(abrupt=ttest_paired([x["abrupt_lick"] for x in r], [x["abrupt_da"] for x in r]),
                change=ttest_paired([x["change_lick"] for x in r], [x["change_da"] for x in r]), per_animal=r)


def stat_T4(M):
    r = [cumsum_stats(m) for m in M.values()]
    return dict(daic_da=ttest1([x["daic_da"] for x in r]), daic_lick=ttest1([x["daic_lick"] for x in r]),
                per_animal=r)


def stat_T6(M):
    r = [cumsum_stats(m) for m in M.values()]
    return dict(abrupt=ttest_paired([x["abrupt_lick"] for x in r], [x["abrupt_da"] for x in r]),
                change=ttest_paired([x["change_lick"] for x in r], [x["change_da"] for x in r]),
                mean_da=ttest1([x["mean_da_ext"] for x in r]), per_animal=r)


# ------------------------------------------------------------------ Test 5 (figS10/omission.m)


def csrw_of(get, sess, subj):
    pav = [(d, n) for d, n, p in sess if p == "Pavlovian"]
    e0 = J.eventfrompavlovian(get(subj, pav[0][1]))
    n = min(len(e0["CS"]), len(e0["fxreward"]))
    return np.unique(e0["CS"][:n][e0["fxreward"][:n] == 1])[0]


def measure_T5(get, sess, subj):
    iM = J.MAIN7.index(subj)
    csrw = csrw_of(get, sess, subj)
    out = {}
    for i, win in enumerate([(3000, 4000), (9000, 10000)]):
        name = [n for d, n, p in sess if d == J.OMIT_DAY[i, iM] and p == "Pavlovian"][0]
        s = get(subj, name)
        e = J.eventfrompavlovian(s)
        cst = e["CStime"][e["CS"] == csrw]
        time, al = J.alignsignal2event(s["t"], s["sig"], cst, [-2000, 12000], 10)
        out["da_" + str(i)] = J.aucsignal(al, time, list(win)) / (win[1] - win[0]) - J.aucsignal(al, time, [-1000, 0])
    return out


def stat_T5(M):
    a = [np.sum(m["da_0"]) for m in M.values()]
    b = [np.sum(m["da_1"]) for m in M.values()]
    return dict(paired=ttest_paired(a, b), old_delay=ttest1([np.mean(m["da_0"]) for m in M.values()]),
                omission=ttest1([np.mean(m["da_1"]) for m in M.values()]))


# ------------------------------------------------------------------ Test 7


def measure_T7(get, sess, subj):
    iM = J.MAIN7.index(subj)
    csrw = csrw_of(get, sess, subj)
    out = {}
    for iC in range(2):
        pav = [(d, n) for d, n, p in sess if p == "Pavlovian" and J.BGD_START[iM, iC] - 2 <= d <= J.BGD_END[iM, iC]]
        aucs, n_pre = [], 0
        for k, (d, n) in enumerate(pav):
            s = get(subj, n)
            e = J.eventfrompavlovian(s)
            a = cs_auc(s, e["CStime"][e["CS"] == csrw])
            aucs.append(a)
            if d < J.BGD_START[iM, iC]:
                n_pre += len(a)
        out["da_%d" % iC] = np.concatenate(aucs)
        out["npre_%d" % iC] = n_pre
    return out


def stat_T7(M):
    v = {0: [], 1: []}
    for m in M.values():
        for iC in range(2):
            data = m["da_%d" % iC][np.arange(140) + m["npre_%d" % iC] - 50]
            v[iC].append(np.cumsum(data)[-1] / np.sum(data[:50]))
    return ttest_paired(v[0], v[1]) | dict(bgd=v[0], ext=v[1])


# ------------------------------------------------------------------ Test 8 (fig5/intermediate_rsp_trialless.m)


def measure_T8(get, sess, subj):
    iM = J.MAIN7.index(subj)
    out = []
    for d, n, p in sess:
        if p != "Poisson" or not (J.POIS_LEARNED[iM] <= d <= J.POIS_END[iM]):
            continue
        s = get(subj, n)
        e = J.eventfrompavlovian(s, trialless=True)
        nn = min(len(e["CS"]), len(e["fxreward"]))
        csrw = np.unique(e["CS"][:nn][e["fxreward"][:nn] == 1])[0]
        nn = min(len(e["CStime"]), len(e["fxrwtime"]))
        cst, fx, first = e["CStime"][:nn], e["fxrwtime"][:nn], e["firstfxlicktime"][:nn]
        delay = np.mean(fx - cst)
        inlick_all = np.r_[first[:-1] < fx[1:], False]
        sel = e["CS"][:nn] == csrw
        tc, cs = J.alignsignal2event(s["t"], s["sig"], cst[sel], [-2000, 10000], 0)
        tr, rw = J.alignsignal2event(s["t"], s["sig"], first[sel], [-2000, 5000], 10)
        base = J.aucsignal(cs, tc, [-500, 0])
        acs = J.aucsignal(cs, tc, [0, 500])
        arw = J.aucsignal(rw, tr, [0, 500], np.isnan(first[sel]))
        cst = cst[sel]
        inlick = inlick_all[sel]
        dif = np.diff(cst)
        inter = np.r_[False, dif < delay]
        ti = np.column_stack([np.r_[inter[1:], False], inter]).astype(int)
        both = np.nonzero(ti.sum(1) == 2)[0]
        ti[both[both - 1 >= 0] - 1, :] = 0
        ti[both[both + 1 < len(ti)] + 1, :] = 0
        ti[ti.sum(1) == 2, :] = 0
        short = np.nonzero(np.r_[False, dif <= 500] | np.r_[False, (dif > 2500) & (dif < 3000)])[0]
        rm = np.r_[short, short - 1]
        ti[rm[rm >= 0], :] = 0
        ti = ti.astype(bool)
        base2 = base.copy()
        base2[ti[:, 1]] = base[ti[:, 0]]
        out.append(dict(da_cs=acs, da_rw=arw, da_base=base2, ti=ti, inlick=inlick))
    return out


def stat_T8(M):
    cue, rew = [], []
    for ses in M.values():
        c = np.full((len(ses), 2), np.nan)
        r = np.full((len(ses), 2), np.nan)
        for i, s in enumerate(ses):
            acs = s["da_cs"] - s["da_base"]
            arw = s["da_rw"] - s["da_base"]
            for ii in range(2):
                c[i, ii] = np.mean(acs[s["ti"][:, ii]]) if s["ti"][:, ii].any() else np.nan
                m = s["ti"][:, ii] & s["inlick"]
                r[i, ii] = np.mean(arw[m]) if m.any() else np.nan
        c, r = np.nanmean(c, 0), np.nanmean(r, 0)
        cue.append(c[1] / c[0]); rew.append(r[1] / r[0])
    return dict(cue_ratio=ttest1(cue, 0.0), reward_ratio=ttest1(rew, 1.0))


# ------------------------------------------------------------------ Test 9 (fig6/backpropagation_pavlovian.m)


def measure_T9(get, sess, subj):
    iM = J.MAIN7.index(subj)
    csrw = csrw_of(get, sess, subj)
    early, late = [], []
    for d, n, p in sess:
        if p != "Pavlovian" or not (J.BP_START[iM] <= d <= J.BP_END[iM]):
            continue
        s = get(subj, n)
        e = J.eventfrompavlovian(s)
        t, al = J.alignsignal2event(s["t"], s["sig"], e["CStime"][e["CS"] == csrw], [-2000, 3000], 0)
        early.append(J.aucsignal(al, t, [0, 1000])); late.append(J.aucsignal(al, t, [2000, 3000]))
    return dict(da_early=np.concatenate(early), da_late=np.concatenate(late))


def stat_T9(M):
    v = [np.mean(m["da_late"][:100] - m["da_early"][:100]) / np.nanmean(m["da_early"][350:400]) for m in M.values()]
    return ttest1(v) | dict(per_animal=[float(x) for x in v])


# ------------------------------------------------------------------ Tests 10, 11 (fig6/*sequential*.py, seconds)


def calc_auc_s(sig, t, ref, win):
    """functions/plot.py calculate_auc (times in the same unit as t)."""
    binsize = np.nanmean(np.diff(t))
    out = []
    for v in ref:
        m = (t >= v + win[0]) & (t <= v + win[1])
        out.append(np.nan if m.sum() == 0 else np.sum(sig[m]) * binsize)
    return np.array(out)


def first_after(ev, ei, target, ref_code, window):
    tt = ev[ei == target]
    out = []
    for v in ev[ei == ref_code]:
        c = tt[(tt > v) & (tt - v <= window)]
        out.append(c[0] if len(c) else np.nan)
    return np.array(out)


def measure_seq(get, sess, subj):
    """Times in ms here; windows converted from the Python scripts' seconds."""
    pav = [(d, n) for d, n, p in sess if p == "Pavlovian"]
    allses = [(d, n) for d, n, p in sess]
    res = dict(da_cs1=[], da_cs2=[], da_cs1_lc=[], da_rw_lc=[])
    for k, (d, n) in enumerate(pav):
        s = get(subj, n)
        et, ei, t, sig = s["eventtime"], s["eventindex"], s["t"], s["sig"]
        cs = et[ei == J.CS1]
        cs1, cs2 = cs[0::2], cs[1::2]
        # backpropagation_sequential.py: AUC 0..1 s, no baseline
        res["da_cs1"].append(calc_auc_s(sig, t, cs1, [0, 1000]))
        res["da_cs2"].append(calc_auc_s(sig, t, cs2, [0, 1000]))
        # learnign_curve_sequential.py: AUC 0..1.5 s minus AUC -1.5..0 s before CS1
        fl = first_after(et, ei, J.LICK, J.FXRW, 5000)
        b = calc_auc_s(sig, t, cs1 - 1500, [0, 1500])
        res["da_cs1_lc"].append(calc_auc_s(sig, t, cs1, [0, 1500]) - b)
        rw = calc_auc_s(sig, t, fl, [0, 1500])
        res["da_rw_lc"].append(rw - b[: len(rw)] if len(rw) <= len(b) else rw[: len(b)] - b)
        if k == 0:
            # measure_inhibition_sequential.py, first Pavlovian session, windows 0..0.5 s
            fl3 = first_after(et, ei, J.LICK, J.FXRW, 3000)
            base = np.nanmean(calc_auc_s(sig, t, cs1 - 500, [0, 500]))
            res["da_inh_rw"] = np.nanmean(calc_auc_s(sig, t, fl3, [0, 500])) - base
            res["da_inh_cs2"] = np.nanmean(calc_auc_s(sig, t, cs2 + 1000, [0, 500])) - base
    return res


def stat_T10(M):
    v = []
    for m in M.values():
        e, l = np.concatenate(m["da_cs1"]), np.concatenate(m["da_cs2"])
        v.append(np.nanmean(l[:50] - e[:50]) / np.nanmean(e[150:200]))
    return ttest1(v) | dict(per_animal=[float(x) for x in v])


def stat_T11(M):
    v = []
    for m in M.values():
        cs1 = [np.nanmean(x) for x in m["da_cs1_lc"]]
        rw = [np.nanmean(x) for x in m["da_rw_lc"]]
        v.append(cs1[-1] / rw[0])
    inh = [m["da_inh_cs2"] / m["da_inh_rw"] for m in M.values()]
    return dict(last_cs1_over_rw=[float(x) for x in v], all_positive=bool(np.all(np.array(v) > 0)),
                inhibition_ratio=[float(x) for x in inh])
