"""Search plan step 3: run the best shared sets on all cohorts (Tests 1 to 11).

Usage: python seqcheck.py CACHE CTX OUT.jsonl IN.jsonl [IN2.jsonl ...] --top 5 --nseeds 2
Takes the best point of each refinement start (or screen point) across inputs,
ranked by tests passed then loss, and scores it with the sequential cohort.
"""
from __future__ import annotations

import argparse
import json

import model_run as R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cache"); ap.add_argument("ctx"); ap.add_argument("out"); ap.add_argument("inp", nargs="+")
    ap.add_argument("--top", type=int, default=5); ap.add_argument("--nseeds", type=int, default=2)
    a = ap.parse_args()
    best = {}
    for f in a.inp:
        for l in open(f):
            r = json.loads(l)
            if "error" in r:
                continue
            key = (f, r.get("start", r["i"]))
            if key not in best or (-r["n_pass"], r["loss"]) < (-best[key]["n_pass"], best[key]["loss"]):
                best[key] = r
    cands = sorted(best.values(), key=lambda r: (-r["n_pass"], r["loss"]))[: a.top]
    ctx = R.Context(a.cache, a.ctx, seq=True)
    for r in cands:
        s = R.score(ctx, r["p"], seq=True, nseeds=a.nseeds)
        rec = dict(src=r.get("start", r["i"]), p=r["p"], n_pass=s["n_pass"], passes=s["passes"], loss=s["loss"],
                   margins=s["margins"], inh=s["inh"], nseeds=a.nseeds)
        with open(a.out, "a") as fh:
            fh.write(json.dumps(rec, default=float) + "\n")
        print(rec["src"], rec["p"].get("cls"), s["n_pass"], {k: v for k, v in s["passes"].items() if not v}, flush=True)


if __name__ == "__main__":
    main()
