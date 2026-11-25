"""Nelder-Mead refinement of the best screened points (METHODS search plan, step 2).

Usage: python refine.py CACHE CTX IN.jsonl OUT.jsonl --cls MS --top 20 --maxfev 60 --workers 3
Categorical parameters stay at the start point's values; continuous ones move in
the unit cube of the search box (clipped). Objective: the screening loss.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp

import numpy as np
from scipy.optimize import minimize

import search as SE


def encode(p, box):
    u = []
    for name, kind, *a in box:
        if kind == "log":
            u.append((np.log(p[name]) - np.log(a[0])) / (np.log(a[1]) - np.log(a[0])))
        elif kind == "lin":
            u.append((p[name] - a[0]) / (a[1] - a[0]))
        else:
            u.append((a[0].index(p[name]) + 0.5) / len(a[0]))
    return np.array(u)


def refine_one(args):
    k, p0, cls, maxfev, out = args
    box = SE.BOXES[cls]
    cont = [i for i, b in enumerate(box) if b[1] != "cat"]
    u0 = encode(p0, box)
    log = []

    def f(z):
        u = u0.copy()
        u[cont] = np.clip(z, 0.0, 1.0)
        p = dict(SE.decode(u, box), cls=cls)
        r = SE.work((k, p))
        log.append(r)
        with open(out, "a") as fh:
            fh.write(json.dumps(dict(r, start=k), default=float) + "\n")
        return r["loss"] - 1000.0 * (r["n_pass"] == len(r["passes"]) and len(r["passes"]) > 0)
    simplex = [u0[cont]] + [np.clip(u0[cont] + 0.08 * np.eye(len(cont))[j], 0, 1) for j in range(len(cont))]
    minimize(f, u0[cont], method="Nelder-Mead",
             options=dict(maxfev=maxfev, initial_simplex=np.array(simplex), xatol=1e-3, fatol=1e-3))
    best = min(log, key=lambda r: (-r["n_pass"], r["loss"]))
    return dict(start=k, best=best)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cache"); ap.add_argument("ctx"); ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--cls", default="MS"); ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--maxfev", type=int, default=60); ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    R = [json.loads(l) for l in open(a.inp)]
    R = [r for r in R if "error" not in r]
    R.sort(key=lambda r: (-r["n_pass"], r["loss"]))
    starts = [(r["i"], r["p"], a.cls, a.maxfev, a.out) for r in R[: a.top]]
    SE.init(a.cache, a.ctx, False)
    with mp.get_context("fork").Pool(a.workers) as pool:
        for res in pool.imap_unordered(refine_one, starts):
            b = res["best"]
            print(res["start"], b["n_pass"], round(b["loss"], 2), flush=True)


if __name__ == "__main__":
    main()
