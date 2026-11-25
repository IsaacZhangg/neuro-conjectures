"""Reproducible check: storage capacity of the Hopfield network.

Claim under test (Amit, Gutfreund & Sompolinsky 1985/1987): a Hopfield network
with Hebbian couplings storing P = alpha*N random binary patterns retrieves them
(stable state close to the stored pattern) only while alpha < alpha_c ~ 0.138,
in the limit N -> infinity, at zero temperature.

What this script does
- Builds Hebbian networks of N binary neurons with P random patterns.
- Starts the network exactly at a stored pattern and runs zero-temperature
  asynchronous dynamics until no neuron wants to flip.
- Records the final overlap m with the starting pattern.
- For each N, estimates the finite-size capacity alpha_c(N) as the load where the
  retrieval probability (m >= 0.9) crosses 1/2, by linear interpolation.
- Extrapolates alpha_c(N) to N -> infinity with a fit linear in N**-0.5.

Run:  python check.py            (full run, writes check-output.json)
      python check.py --quick    (smoke test for CI; does not set the status)
Requires Python 3.10+ and NumPy. Fully deterministic given the seed.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time

import numpy as np

PREDICTED_ALPHA_C = 0.138  # AGS replica-symmetric value, 0.138 (often quoted 0.1379)
TOLERANCE = 0.02           # pass if the extrapolated capacity is within +/- 0.02


def run_trial(N: int, P: int, starts: int, rng: np.random.Generator) -> list[float]:
    xi = rng.choice(np.array([-1, 1], dtype=np.int8), size=(P, N)).astype(np.float64)
    W = (xi.T @ xi) / N
    np.fill_diagonal(W, 0.0)
    overlaps = []
    for mu in range(min(starts, P)):
        s = xi[mu].copy()
        h = W @ s
        # Asynchronous zero-temperature dynamics: repeatedly flip one randomly
        # chosen unstable neuron. Energy strictly decreases, so this terminates.
        for _ in range(50 * N):
            unstable = np.flatnonzero(s * h < 0)
            if unstable.size == 0:
                break
            i = rng.choice(unstable)
            s[i] = -s[i]
            h += 2.0 * s[i] * W[:, i]
        else:
            raise RuntimeError("dynamics did not converge")
        overlaps.append(float(s @ xi[mu]) / N)
    return overlaps


def crossing(alphas: list[float], probs: list[float], level: float = 0.5) -> float | None:
    for (a0, p0), (a1, p1) in zip(zip(alphas, probs), zip(alphas[1:], probs[1:])):
        if p0 >= level > p1:
            return a0 + (p0 - level) * (a1 - a0) / (p0 - p1)
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--seed", type=int, default=20261125)
    args = parser.parse_args()

    if args.quick:
        sizes, networks, starts = [250, 500], 4, 5
    else:
        sizes, networks, starts = [500, 1000, 2000, 4000], 12, 10
    alphas = [round(0.10 + 0.01 * k, 2) for k in range(15)]  # 0.10 .. 0.24

    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    table = {}
    for N in sizes:
        rows = []
        for a in alphas:
            P = max(1, round(a * N))
            ms = []
            for _ in range(networks):
                ms.extend(run_trial(N, P, starts, rng))
            ms = np.array(ms)
            rows.append({
                "alpha": a,
                "P": P,
                "samples": int(ms.size),
                "mean_overlap": round(float(ms.mean()), 4),
                "retrieval_probability": round(float((ms >= 0.9).mean()), 4),
            })
        ac = crossing([r["alpha"] for r in rows], [r["retrieval_probability"] for r in rows])
        table[N] = {"rows": rows, "alpha_c_finite_N": None if ac is None else round(ac, 4)}
        print(f"N={N:5d}  alpha_c(N) ~ {table[N]['alpha_c_finite_N']}", flush=True)

    pts = [(N, v["alpha_c_finite_N"]) for N, v in table.items() if v["alpha_c_finite_N"] is not None]
    extrapolated = None
    if len(pts) >= 2:
        x = np.array([N ** -0.5 for N, _ in pts])
        y = np.array([a for _, a in pts])
        slope, intercept = np.polyfit(x, y, 1)
        extrapolated = round(float(intercept), 4)

    passed = extrapolated is not None and abs(extrapolated - PREDICTED_ALPHA_C) <= TOLERANCE
    if args.quick:
        # Quick mode is a smoke test for CI: it only checks that the pipeline runs
        # and finds a transition. The status comes from the full run.
        ok = all(v["alpha_c_finite_N"] is not None for v in table.values())
        print(f"quick smoke test: transition found for every N = {ok} (extrapolation {extrapolated}, not used for status)")
        return 0 if ok else 1
    out = {
        "result": "hopfield-capacity",
        "mode": "quick" if args.quick else "full",
        "seed": args.seed,
        "predicted_alpha_c": PREDICTED_ALPHA_C,
        "tolerance": TOLERANCE,
        "retrieval_threshold_overlap": 0.9,
        "finite_size": {str(N): v for N, v in table.items()},
        "extrapolation": "alpha_c(N) = alpha_c + c * N**-0.5, least squares",
        "alpha_c_extrapolated": extrapolated,
        "passed": bool(passed),
        "runtime_seconds": round(time.time() - t0, 1),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
    }
    print(f"extrapolated alpha_c = {extrapolated}  (predicted {PREDICTED_ALPHA_C} +/- {TOLERANCE})  passed={passed}")
    if not args.quick:
        with open("check-output.json", "w") as f:
            json.dump(out, f, indent=2)
            f.write("\n")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
