"""Summarise a search JSONL: pass rate per test, joint passes, best points."""
import json
import sys
import collections
import itertools

import numpy as np

R = [json.loads(l) for l in open(sys.argv[1])]
ok = [r for r in R if "error" not in r]
print(len(R), "points,", len(R) - len(ok), "diverged or failed")
tests = sorted({k for r in ok for k in r["passes"]}, key=lambda s: int(s[1:]))
print("pass rate:", {t: round(np.mean([r["passes"].get(t, False) for r in ok]), 3) for t in tests})
print("n_pass histogram:", dict(sorted(collections.Counter(r["n_pass"] for r in ok).items())))
print("never jointly passed (pairs):", [a + "+" + b for a, b in itertools.combinations(tests, 2)
                                         if not any(r["passes"].get(a) and r["passes"].get(b) for r in ok)])
best = sorted(ok, key=lambda r: (-r["n_pass"], r["loss"]))[: int(sys.argv[2]) if len(sys.argv) > 2 else 5]
for r in best:
    print(r["n_pass"], round(r["loss"], 1), "fails", [t for t in tests if not r["passes"].get(t)],
          {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r["p"].items()})
# best margin reached per criterion
mk = sorted({k for r in ok for k in r.get("margins", {})})
print("best margin per criterion:", {k: round(max(r["margins"].get(k, -1e9) for r in ok), 2) for k in mk})
