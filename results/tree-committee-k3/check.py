"""Certified check of the K=3 tree committee machine lower bound (Theorem 1 in PROOF.md).

Every inequality the proof needs from numerics is evaluated in Arb ball
arithmetic (python-flint), so each printed comparison is rigorous. Exit 0 if
all of them hold, 1 otherwise. Writes check-output.json.

    python check.py           # full run (same as quick; it takes about a second)
    python check.py --quick
"""
from __future__ import annotations

import json
import platform
import sys
import time
from fractions import Fraction as F
from pathlib import Path

import flint
from flint import arb, ctx

ctx.prec = 200
SQ2 = arb(2).sqrt()
SQ2PI = (2 * arb.pi()).sqrt()


def q(v) -> arb:
    v = F(v)
    return arb(v.numerator) / arb(v.denominator)


def Phi(t):  # standard normal CDF
    return (-t / SQ2).erfc() / 2


def phi(t):  # standard normal density
    return (-(t * t) / 2).exp() / SQ2PI


def I2(lo, hi, c):  # int_lo^hi (c-g)^2 phi(g) dg
    P = Phi(hi) - Phi(lo)
    M = phi(lo) - phi(hi)
    S = P - (hi * phi(hi) - lo * phi(lo))
    return c * c * P - 2 * c * M + S


def I2_from_minus_inf(hi, c):  # int_-inf^hi (c-g)^2 phi(g) dg
    P = Phi(hi)
    S = P - hi * phi(hi)
    return c * c * P + 2 * c * phi(hi) + S


def soft(x, rho, k):
    """Soft population, displacement D(g) = clamp(k-g, x*rho, x).
    Returns E[D^2], E[(g+D)D] and P(g + D(g) < k) = Phi(k-x)."""
    A, B = k - x, k - x * rho
    c = x * x * Phi(A) + I2(A, B, k) + (x * rho) ** 2 * (1 - Phi(B))
    return c, c - (Phi(B) - Phi(A)), Phi(A)


def hard(x, rho, k):
    """Hard population, displacement D(g) = max(k-g, x*rho). Returns E[D^2], E[(g+D)D]."""
    B = k - x * rho
    c = I2_from_minus_inf(B, k) + (x * rho) ** 2 * (1 - Phi(B))
    return c, c - Phi(B)


def stage(LH, LR, rho, k, width=F(1, 10**6)):
    """Bracket the root x* of LH*E[D_H^2] + LR*E[D_R^2] = 1 (the cost is increasing
    in x), and enclose the ball-multiplier quantity P on the whole bracket."""
    def cost(x):
        c = LR * soft(x, rho, k)[0]
        if LH != 0:
            c += LH * hard(x, rho, k)[0]
        return c

    lo, hi = F(1, 100), F(5)
    assert bool(cost(q(lo)) < 1) and bool(cost(q(hi)) > 1)
    while hi - lo > width:
        mid = (lo + hi) / 2
        v = cost(q(mid)) - 1
        if v > 0:
            hi = mid
        elif v < 0:
            lo = mid
        else:
            raise RuntimeError("cannot decide sign; raise precision")
    ball = arb((q(lo) + q(hi)) / 2, (q(hi) - q(lo)) / 2)
    P = LR * soft(ball, rho, k)[1]
    if LH != 0:
        P += LH * hard(ball, rho, k)[1]
    err_upper = soft(q(lo), rho, k)[2]  # Phi(k - x) is decreasing in x
    return dict(x_lo=lo, x_hi=hi, cost_lo=cost(q(lo)), cost_hi=cost(q(hi)), P=P, err_upper=err_upper)


# Parameters of Theorem 1 (all exact rationals)
ALPHA = F(31, 20)      # 1.55
KAPPA = F(1, 1000)     # margin used for "correct"
RHO1 = F(59, 250)      # stage 1 loss ratio
RHO2 = F(7, 25)        # stage 2 loss ratio (0.28), used for both populations
A_UP = F(975253, 10**6)  # rational upper bound for the stage-1 error load


def main() -> int:
    t0 = time.time()
    beta = 3 * q(ALPHA)
    k = q(KAPPA)
    s1 = stage(0, beta, q(RHO1), k)
    load1 = beta * s1["err_upper"]
    a = q(A_UP)
    r = beta - a
    s2 = stage(a, r, q(RHO2), k)
    final = a + r * s2["err_upper"]
    checks = {
        "stage1: P > 0 on the root bracket": bool(s1["P"] > 0),
        "stage1: beta*e1_upper < A_UP": bool(load1 < a),
        "stage2: A_UP < 1 (needed for the multiplier bound)": bool(a < 1),
        "stage2: P > 0 on the root bracket": bool(s2["P"] > 0),
        "stage3: A_UP + (beta-A_UP)*e2_upper < 2 (Cover)": bool(final < 2),
    }
    ok = all(checks.values())
    s = lambda v: v.str(15, radius=True) if isinstance(v, arb) else str(v)
    out = {
        "claim": "K=3 tree committee machine, Gaussian patterns: all m = floor(alpha n) patterns storable w.h.p. for alpha = 31/20",
        "parameters": {"alpha": str(ALPHA), "kappa": str(KAPPA), "rho1": str(RHO1), "rho2": str(RHO2), "a_up": str(A_UP)},
        "stage1": {kk: s(v) for kk, v in s1.items()} | {"beta*e1_upper": s(load1)},
        "stage2": {kk: s(v) for kk, v in s2.items()},
        "stage3_load": s(final),
        "checks": checks,
        "pass": ok,
        "mode": "quick" if "--quick" in sys.argv else "full",
        "software": {"python": platform.python_version(), "python-flint": flint.__version__, "arb_prec_bits": ctx.prec},
        "runtime_s": round(time.time() - t0, 3),
    }
    for kk, v in checks.items():
        print(("PASS " if v else "FAIL ") + kk)
    print("stage-3 load", s(final), "(must be < 2)")
    Path(__file__).with_name("check-output.json").write_text(json.dumps(out, indent=2) + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
