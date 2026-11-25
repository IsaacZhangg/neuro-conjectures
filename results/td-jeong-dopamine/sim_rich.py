"""TD-BS and TD-SR simulators (METHODS.md, Addendum A).

Same event streams, delta rule, eligibility trace, kernel and scaling as TD-MS
(sim.py). Implementation choices the addendum leaves open are listed in
METHODS.md, deviation 7.
"""
from __future__ import annotations

import numpy as np
from numba import njit

import sim as S

K = 24  # wait sub-states per cue
WBIN = 5  # bins of 0.1 s per sub-state (0.5 s)
NC = 2  # cue types: A (stream type 0) and B (stream type 2)


@njit(cache=True)
def run_session_bs(evb, evt, rw, laser, nb, bias, w, gam, lam, alpha, eta, inh, kappa, cN, cR, cH, out):
    nfeat = w.shape[0]
    x = np.zeros(nfeat)
    xp = np.zeros(nfeat)
    e = np.zeros(nfeat)
    on = np.full(64, -1, np.int64)  # open instances: onset bin
    ty = np.zeros(64, np.int64)
    nxt = 0
    k = 0
    vprev = 0.0
    nev = evb.shape[0]
    for b in range(nb):
        while k < nev and evb[k] <= b:
            t = evt[k]
            if t == 0 or t == 2:
                on[nxt % 64] = evb[k]
                ty[nxt % 64] = t // 2
                nxt += 1
            k += 1
        # rewards close the instance that best predicts a reward now
        nr = int(rw[b] + 0.5)
        for _ in range(nr):
            best, bs, recent = -1, 0.0, -1
            for j in range(64):
                if on[j] < 0:
                    continue
                c = ty[j]
                kk = (b - on[j]) // WBIN
                if kk >= K:
                    continue
                if recent < 0 or on[j] >= on[recent]:
                    recent = j
                if cR[c] > 1e-12:
                    sc = (cR[c] / cN[c]) * cH[c, kk] / cR[c]
                    if sc > bs:
                        bs, best = sc, j
            if best < 0:
                best = recent
            if best >= 0:
                c = ty[best]
                kk = (b - on[best]) // WBIN
                cN[c] = (1 - kappa) * cN[c] + 1.0
                cR[c] = (1 - kappa) * cR[c] + 1.0
                for i in range(K):
                    cH[c, i] *= (1 - kappa)
                cH[c, kk] += 1.0
                on[best] = -1
        # time-outs close unrewarded instances
        for j in range(64):
            if on[j] >= 0 and (b - on[j]) // WBIN >= K:
                c = ty[j]
                cN[c] = (1 - kappa) * cN[c] + 1.0
                cR[c] = (1 - kappa) * cR[c]
                for i in range(K):
                    cH[c, i] *= (1 - kappa)
                on[j] = -1
        for i in range(nfeat):
            x[i] = 0.0
        tot = 0.0
        for j in range(64):
            if on[j] < 0:
                continue
            c = ty[j]
            kk = (b - on[j]) // WBIN
            pi = cR[c] / cN[c]
            surv = 1.0
            if cR[c] > 1e-12:
                acc = 0.0
                for i in range(kk):
                    acc += cH[c, i]
                surv = max(0.0, 1.0 - acc / cR[c])
            num = pi * surv
            den = num + 1.0 - pi
            pj = num / den if den > 1e-12 else 0.0
            x[c * K + kk] += pj
            tot += pj
        x[NC * K] = min(1.0, max(0.0, 1.0 - tot))
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
def _phi(b, occ, table, use_off, phi):
    nf_t = table.shape[1]
    for i in range(phi.shape[0]):
        phi[i] = 0.0
    for t in range(S.NT):
        for j in range(64):
            s0 = occ[t, j]
            if s0 < 0:
                continue
            a = b - s0
            if a >= S.MAXAGE:
                occ[t, j] = -1
                continue
            base = t * nf_t
            for i in range(nf_t):
                phi[base + i] += table[a, i]


@njit(cache=True)
def run_session_sr(evb, evt, rw, laser, nb, table, use_off, bias, w, M, gam, lam, alpha, alpha_m, eta, inh, out):
    nf = M.shape[0]
    phi = np.zeros(nf)
    phip = np.zeros(nf)
    psi = np.zeros(nf)
    psip = np.zeros(nf)
    e = np.zeros(nf)
    err = np.zeros(nf)
    occ = np.full((S.NT, 64), -1, np.int64)
    nocc = np.zeros(S.NT, np.int64)
    nzp = np.zeros(nf, np.int64)
    nz = np.zeros(nf, np.int64)
    k = 0
    vprev = 0.0
    nev = evb.shape[0]
    nnzp = 0
    for b in range(nb):
        while k < nev and evb[k] <= b:
            t = evt[k]
            if use_off or (t != 1 and t != 3):
                occ[t, nocc[t] % 64] = evb[k]
                nocc[t] += 1
            k += 1
        _phi(b, occ, table, use_off, phi)
        if bias:
            phi[nf - 1] = 1.0
        nnz = 0
        for i in range(nf):
            if abs(phi[i]) > 1e-7:
                nz[nnz] = i
                nnz += 1
        # psi = M phi
        for r in range(nf):
            acc = 0.0
            for q in range(nnz):
                acc += M[r, nz[q]] * phi[nz[q]]
            psi[r] = acc
        v = 0.0
        for i in range(nf):
            v += w[i] * psi[i]
        if b == 0:
            vprev = v
        d = rw[b] + gam * v - vprev - inh * laser[b]
        out[b] = d
        g = d if d >= 0 else eta * d
        for i in range(nf):
            e[i] = gam * lam * e[i] + psip[i]
            w[i] += alpha * g * e[i]
        # successor update: M += a_M (phi_{t-1} + gam M phi_t - M phi_{t-1}) phi_{t-1}^T
        if alpha_m > 0 and nnzp > 0:
            for r in range(nf):
                err[r] = phip[r] + gam * psi[r] - psip[r]
            for q in range(nnzp):
                c = nzp[q]
                f = alpha_m * phip[c]
                for r in range(nf):
                    M[r, c] += f * err[r]
        for i in range(nf):
            phip[i] = phi[i]
            psip[i] = psi[i]
        for q in range(nnz):
            nzp[q] = nz[q]
        nnzp = nnz
        vprev = v


def run_bs(animal, p, inh=0.0, state=None):
    nfeat = NC * K + 1 + (1 if p["bias"] else 0)
    if state is None:
        state = dict(w=np.full(nfeat, p["w0"]), N=np.ones(NC), R=np.zeros(NC), H=np.zeros((NC, K)))
    gam = p["gamma_s"] ** S.DT
    out = {}
    for d, n, st in animal.streams:
        delta = np.zeros(st["nb"])
        run_session_bs(st["evb"], st["evt"], st["rw"], st["laser"], st["nb"], bool(p["bias"]), state["w"], gam,
                       p["lam"], p["alpha"], p["eta"], inh, p["kappa"], state["N"], state["R"], state["H"], delta)
        if n is not None:
            out[n] = (st, delta)
    return out, state


def run_sr(animal, p, table, inh=0.0, state=None):
    nf = S.NT * table.shape[1] + (1 if p["bias"] else 0)
    if state is None:
        state = dict(w=np.full(nf, p["w0"]), M=np.eye(nf))
    gam = p["gamma_s"] ** S.DT
    out = {}
    for d, n, st in animal.streams:
        delta = np.zeros(st["nb"])
        run_session_sr(st["evb"], st["evt"], st["rw"], st["laser"], st["nb"], table, bool(p["offset"]),
                       bool(p["bias"]), state["w"], state["M"], gam, p["lam"], p["alpha"], p["alpha_m"], p["eta"],
                       inh, delta)
        if n is not None:
            out[n] = (st, delta)
    return out, state
