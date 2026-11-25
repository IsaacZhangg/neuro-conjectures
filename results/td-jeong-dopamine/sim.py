"""TD-MS simulator (METHODS.md, "Model family TD-MS").

Each animal's sessions are turned into event streams on a 0.1 s grid. The model
runs through all sessions in order, weights carried over, traces reset at each
session start. Output per session: TD error per bin, and simulated dF/F at the
data's own timestamps (kernel h, scale 1; the scale is applied later).
"""
from __future__ import annotations

import numpy as np
from numba import njit

import jeong as J

DT = 0.1  # s
NT = 5  # event types: 0 A on, 1 A off, 2 B on, 3 B off, 4 reward
MAXAGE = 600  # 60 s superposition cut
KLEN = 40  # kernel taps (4 s)
TAU_D, TAU_R = 0.5, 0.05


def kernel_value(tau, tau_d=TAU_D, tau_r=TAU_R):
    tp = np.log(tau_d / tau_r) * tau_d * tau_r / (tau_d - tau_r)
    pk = np.exp(-tp / tau_d) - np.exp(-tp / tau_r)
    v = (np.exp(-tau / tau_d) - np.exp(-tau / tau_r)) / pk
    return np.where(tau >= 0, v, 0.0)


# ------------------------------------------------------------------ session -> event stream
def prepare(s: dict, csA: int | None, kind: str, cs_dur_ms: float | None):
    """kind: 'rr', 'pav', 'pois', 'seq'. csA: code of sound A (CS+), or None."""
    et = s["eventtime"] / 1000.0
    ei = s["eventindex"]
    ns = s["nonsolenoidflag"]
    tt = s["t"] / 1000.0
    good = ~np.isnan(tt)
    t0 = min(np.nanmin(tt) if good.any() else 0.0, et.min()) - 1.0
    t1 = max(np.nanmax(tt) if good.any() else 0.0, et.max()) + 5.0
    nb = int(np.ceil((t1 - t0) / DT)) + 1
    evb, evt = [], []
    lick = et[ei == J.LICK]

    def add(times, typ):
        for x in np.atleast_1d(times):
            if not np.isnan(x):
                evb.append(int((x - t0) / DT)); evt.append(typ)

    # rewards: first lick within 5 s after a delivered solenoid, else the solenoid time
    sol = np.r_[et[(ei == J.FXRW) & (ns == 0)], et[ei == J.BGD]]
    rw = np.zeros(nb, np.float64)
    for x in sol:
        c = lick[(lick > x) & (lick <= x + 5.0)]
        tr = c[0] if len(c) else x
        b = int((tr - t0) / DT)
        rw[b] += 1.0
        add(tr, 4)
    if kind == "seq":
        cs = et[ei == J.CS1]
        add(cs[0::2], 0); add(cs[1::2], 2)
    else:
        for code, on in ((csA, 0), ({J.CS1: J.CS2, J.CS2: J.CS1}.get(csA), 2)):
            if code is None:
                continue
            x = et[ei == code]
            add(x, on)
            if cs_dur_ms is not None:
                add(x + cs_dur_ms / 1000.0, on + 1)
    laser = np.zeros(nb, np.float64)
    if kind == "seq":
        lz = et[ei == 31]
        fx = et[ei == J.FXRW]
        for x in lz:  # laser from CS2 until reward (paper); end = next cued solenoid
            nxt = fx[fx >= x]
            end = nxt[0] if len(nxt) else x + 1.5
            laser[int((x - t0) / DT):int((end - t0) / DT)] = 1.0
    order = np.argsort(np.array(evb, dtype=np.int64), kind="stable")
    return dict(nb=nb, t0=t0, evb=np.array(evb, np.int64)[order], evt=np.array(evt, np.int64)[order],
                rw=rw, laser=laser, tt=tt)


def ms_table(m, sigma, d_per_s, csc=False):
    """f[age, i] for age in bins 0..MAXAGE-1."""
    age = np.arange(MAXAGE) * DT
    if csc:
        nf = 100
        f = np.zeros((MAXAGE, nf))
        for a in range(min(nf, MAXAGE)):
            f[a, a] = 1.0
        return f
    y = d_per_s ** age
    mu = np.arange(1, m + 1) / m
    return y[:, None] * np.exp(-(y[:, None] - mu[None, :]) ** 2 / (2 * sigma ** 2))


@njit(cache=True)
def run_session(evb, evt, rw, laser, nb, table, use_off, bias, w, gam, lam, alpha, eta, inh, out):
    nf_t = table.shape[1]
    nfeat = w.shape[0]
    x = np.zeros(nfeat)
    xp = np.zeros(nfeat)
    e = np.zeros(nfeat)
    # active occurrences: start bins per type (ring of 64)
    occ = np.full((NT, 64), -1, np.int64)
    nocc = np.zeros(NT, np.int64)
    k = 0
    vprev = 0.0
    nev = evb.shape[0]
    for b in range(nb):
        while k < nev and evb[k] <= b:
            ty = evt[k]
            if use_off or (ty != 1 and ty != 3):
                occ[ty, nocc[ty] % 64] = evb[k]
                nocc[ty] += 1
            k += 1
        for i in range(nfeat):
            x[i] = 0.0
        for ty in range(NT):
            for j in range(64):
                s0 = occ[ty, j]
                if s0 < 0:
                    continue
                a = b - s0
                if a >= MAXAGE:
                    occ[ty, j] = -1
                    continue
                base = ty * nf_t
                for i in range(nf_t):
                    x[base + i] += table[a, i]
        if bias:
            x[nfeat - 1] = 1.0
        v = 0.0
        for i in range(nfeat):
            v += w[i] * x[i]
        if b == 0:
            vprev = v
        d = rw[b] + gam * v - vprev - inh * laser[b]
        out[b] = d
        g = d if d >= 0 else eta * d
        for i in range(nfeat):
            e[i] = gam * lam * e[i] + xp[i]
            w[i] += alpha * g * e[i]
            xp[i] = x[i]
        vprev = v


@njit(cache=True)
def trace_at(delta, t0, dt, tt, tau_d, tau_r, norm):
    """Exact y(t) = sum_b delta_b h(t - t_b), h = (exp(-tau/tau_d) - exp(-tau/tau_r)) / norm,
    delta_b placed at t_b = t0 + b dt, via two recursive exponential filters."""
    nb = delta.shape[0]
    A = np.zeros(nb)
    B = np.zeros(nb)
    ad, ar = np.exp(-dt / tau_d), np.exp(-dt / tau_r)
    a = 0.0
    c = 0.0
    for b in range(nb):
        a = a * ad + delta[b]
        c = c * ar + delta[b]
        A[b] = a
        B[b] = c
    n = tt.shape[0]
    y = np.full(n, np.nan)
    for j in range(n):
        t = tt[j]
        if np.isnan(t):
            continue
        b = int((t - t0) / dt)
        if b < 0:
            y[j] = 0.0
            continue
        if b >= nb:
            b = nb - 1
        tau = t - (t0 + b * dt)
        y[j] = (A[b] * np.exp(-tau / tau_d) - B[b] * np.exp(-tau / tau_r)) / norm
    return y


def _kernel_norm(tau_d=TAU_D, tau_r=TAU_R):
    tp = np.log(tau_d / tau_r) * tau_d * tau_r / (tau_d - tau_r)
    return float(np.exp(-tp / tau_d) - np.exp(-tp / tau_r))


KNORM = _kernel_norm()


class Animal:
    """Event streams for one animal's whole history (data-only, model-agnostic)."""

    def __init__(self, cache, subj, fill_gaps=False, tau_d=TAU_D):
        self.subj = subj
        self.tau_d = tau_d
        self.sess = J.list_sessions(cache, subj)
        self.raw = {n: J.load(cache, subj, n) for _, n, _ in self.sess}
        csA = None
        for d, n, p in self.sess:
            if p == "Pavlovian" and "stGtACR" not in subj:
                e = J.eventfrompavlovian(self.raw[n])
                nn = min(len(e["CS"]), len(e["fxreward"]))
                cs = np.unique(e["CS"][:nn][e["fxreward"][:nn] == 1])[0]
                csA = J.CS1 if cs == 1 else J.CS2
                break
        self.csA = csA
        self.streams = []
        for d, n, p in self.sess:
            s = self.raw[n]
            if p == "RandomRewards":
                kind, dur = "rr", None
            elif "stGtACR" in subj:
                kind, dur = "seq", None
            elif p == "Poisson":
                kind, dur = "pois", 250.0
            else:
                kind = "pav"
                dur = J.eventfrompavlovian(s)["CSduration"] if np.any(s["eventindex"] == J.FXRW) else None
            self.streams.append((d, n, prepare(s, csA, kind, dur)))
        if fill_gaps:
            self.streams = fill_short_gaps(self.streams, self.sess)

    def run(self, p, table, inh=0.0, w=None):
        nf_t = table.shape[1]
        nfeat = NT * nf_t + (1 if p["bias"] else 0)
        if w is None:
            w = np.full(nfeat, p["w0"])
        gam = p["gamma_s"] ** DT
        out = {}
        for d, n, st in self.streams:
            delta = np.zeros(st["nb"])
            run_session(st["evb"], st["evt"], st["rw"], st["laser"], st["nb"], table, bool(p["offset"]),
                        bool(p["bias"]), w, gam, p["lam"], p["alpha"], p["eta"], inh, delta)
            if n is not None:
                out[n] = (st, delta)
        return out, w

    def model_sessions(self, out):
        """Session dicts with the model's trace in place of dF/F (same timestamps, NaN mask)."""
        res = {}
        for n, (st, delta) in out.items():
            s = self.raw[n]
            y = trace_at(delta, st["t0"], DT, st["tt"], self.tau_d, TAU_R, _kernel_norm(self.tau_d))
            y[np.isnan(s["sig"])] = np.nan
            res[n] = dict(s, sig=y)
        return res


def fill_short_gaps(streams, sess, maxgap=5):
    """Robustness rule: gaps of <= maxgap days inside one protocol block get copies
    of the preceding session (n=None marks a copy, not analysed)."""
    out = []
    prot = {n: p for _, n, p in sess}
    for i, (d, n, st) in enumerate(streams):
        out.append((d, n, st))
        if i + 1 < len(streams):
            d2, n2, _ = streams[i + 1]
            gap = d2 - d - 1
            if 0 < gap <= maxgap and prot[n] == prot[n2]:
                for _ in range(gap):
                    out.append((d, None, st))
    return out
