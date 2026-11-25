"""Sobol screen of TD-MS over the METHODS box (alpha box per deviation 6).

Usage: python search.py CACHE CTX OUT.jsonl [--n 4096] [--workers 3] [--seq]
Screening: Tests 1-9 (Tests 10-11 for finalists), 2 noise seeds, 6 Weibull starts.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import warnings

import numpy as np
from scipy.stats import qmc

BOX_MS = [("alpha", "log", 1e-5, 0.3), ("gamma_s", "lin", 0.5, 0.999), ("lam", "lin", 0.0, 0.99),
       ("m", "cat", [5, 10, 20]), ("sigma", "lin", 0.02, 0.2), ("d", "lin", 0.3, 0.99),
       ("w0", "lin", 0.0, 0.5), ("bias", "cat", [0, 1]), ("offset", "cat", [0, 1]),
       ("eta", "lin", 0.02, 1.0), ("csc", "cat", [False, True])]
BOX_BS = [("alpha", "log", 1e-5, 0.3), ("gamma_s", "lin", 0.5, 0.999), ("lam", "lin", 0.0, 0.99),
          ("eta", "lin", 0.02, 1.0), ("w0", "lin", 0.0, 0.5), ("kappa", "log", 0.01, 0.5), ("bias", "cat", [0, 1])]
BOX_SR = BOX_MS + [("alpha_m", "log", 1e-4, 0.1)]
BOXES = dict(MS=BOX_MS, BS=BOX_BS, SR=BOX_SR)
BOX = BOX_MS


def decode(u, box=None):
    box = box or BOX
    p = {}
    for x, (name, kind, *a) in zip(u, box):
        if kind == "log":
            p[name] = float(np.exp(np.log(a[0]) + x * (np.log(a[1]) - np.log(a[0]))))
        elif kind == "lin":
            p[name] = float(a[0] + x * (a[1] - a[0]))
        else:
            ch = a[0]
            p[name] = ch[min(int(x * len(ch)), len(ch) - 1)]
    return p


_ctx = None


def init(cache, ctxfile, seq):
    global _ctx
    warnings.filterwarnings("ignore")
    import tests as T
    T.FIT_STARTS = 6
    import model_run as R
    _ctx = R.Context(cache, ctxfile, seq=seq)  # stats come from ctxfile; only needed animals load
    import gc
    for s, a in _ctx.animals.items():  # free what the measures never read
        a.raw = {n: v for n, v in a.raw.items() if (s, n) in _ctx.used}
        a.streams = [(d, n, dict(st, tt=st["tt"] if (s, n) in _ctx.used else None)) for d, n, st in a.streams]
    gc.collect()
    _ctx.seq = seq


def work(item):
    i, p = item
    import model_run as R
    try:
        r = R.score(_ctx, p, seq=_ctx.seq, nseeds=2)
        return dict(i=i, p=p, n_pass=r["n_pass"], passes=r["passes"], loss=r["loss"], margins=r["margins"],
                    sec=r["seconds"])
    except Exception as ex:  # divergence or degenerate output
        return dict(i=i, p=p, n_pass=0, passes={}, loss=1e9, error=repr(ex)[:200])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cache"); ap.add_argument("ctx"); ap.add_argument("out")
    ap.add_argument("--n", type=int, default=4096); ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--seq", action="store_true"); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cls", default="MS", choices=["MS", "BS", "SR"])
    a = ap.parse_args()
    box = BOXES[a.cls]
    U = qmc.Sobol(d=len(box), scramble=True, seed=a.seed).random(a.n)
    done = set()
    if os.path.exists(a.out):
        done = {json.loads(l)["i"] for l in open(a.out)}
    items = [(i, dict(decode(u, box), cls=a.cls)) for i, u in enumerate(U) if i not in done]
    init(a.cache, a.ctx, a.seq)  # load once; forked workers share the arrays copy-on-write
    with mp.get_context("fork").Pool(a.workers) as pool, open(a.out, "a") as f:
        for r in pool.imap_unordered(work, items, chunksize=1):
            f.write(json.dumps(r, default=float) + "\n"); f.flush()


if __name__ == "__main__":
    main()
