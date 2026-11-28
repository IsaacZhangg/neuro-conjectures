"""Check for results/dense-am-capacity.

Part 1 (decides the status together with PROOF.md): certify, with Arb interval
arithmetic, every numerical constant the proof uses.
Part 2 (cross-check): recompute LM's union-bound threshold directly from its
definition on grids (REM free energy by grid maximization, no closed forms)
and compare with the closed form min(alpha_a, log2/4) of PROOF.md Lemma 0.
Part 3 (support only, full mode): brute-force P(all retrieved) at
lambda = inf for small N. Never used to decide anything.

Usage: python check.py [--quick]
Exits 0 if every certified check and the cross-check pass, else 1.
The full run writes check-output.json; quick mode writes nothing.
"""
import json
import math
import platform
import sys
import time

import numpy as np
import flint
from flint import arb, ctx

ctx.prec = 256
QUICK = "--quick" in sys.argv
SEED = 20261125
t_start = time.time()
results = {"certified": {}, "crosscheck": {}, "simulation": {}}
ok = True


def record(name, passed, **info):
    global ok
    ok = ok and bool(passed)
    results["certified"][name] = {"pass": bool(passed), **{k: str(v) for k, v in info.items()}}
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: " + ", ".join(f"{k}={v}" for k, v in info.items()))


# ---------- Part 1: certified constants ----------
L2 = arb(2).log()
alpha_b = L2 / 4


def I(s):
    return (s - 1 - s.log()) / 2


I_half = I(arb(1) / 2)
record("I(1/2) = (log2 - 1/2)/2 > 0", I_half > 0, value=I_half)
# I(s) + s/2 is convex with derivative 1 - 1/(2s): minimum at s = 1/2.
record("min_s I(s)+s/2 = log2/2 at s=1/2",
       abs(I_half + arb(1) / 4 - L2 / 2) < arb(10) ** -60, value=I_half + arb(1) / 4)


def rho_star(r):
    return (-(1 + r) + ((1 + r) ** 2 + 4 * r).sqrt()) / 2


def g(rho, r):
    return (1 - rho * rho).log() / 2 - r / (1 + rho)


r = arb(1) / 2
rs = rho_star(r)
# g(., r) is concave on (-1, 1) (sum of two concave functions) and rs solves
# g'(rho) = -rho/(1-rho^2) + r/(1+rho)^2 = 0, so h(r) = sup g = g(rs, r).
gp = -rs / (1 - rs * rs) + r / (1 + rs) ** 2
record("rho* is the stationary point of g (r=1/2)", abs(gp) < arb(10) ** -60, rho_star=rs)
h_half = g(rs, r)
alpha_star = I_half - h_half - arb(1) / 4
record("alpha_* = I(1/2) - h(1/2) - 1/4 > alpha_b (Section 2 window non-empty)",
       alpha_star > alpha_b, h_half=h_half, alpha_star=alpha_star, alpha_b=alpha_b)
record("alpha_b - I(1/2) > 0 (shell count grows)", alpha_b - I_half > 0, value=alpha_b - I_half)
inv2e = 1 / (2 * arb(1).exp())
record("alpha_a = 1/(2e) > alpha_b for lambda >= 1", inv2e > alpha_b, value=inv2e)


def alpha_a_arb(lam):
    lam = arb(lam)
    a = (1 - lam) ** 2
    c = lam * (1 - lam / 2)
    w = (-a / arb(1).exp()).lambertw()  # principal branch W0
    return -c * w / a


lo, hi = alpha_a_arb("0.700"), alpha_a_arb("0.701")
record("lambda_x in (0.700, 0.701): alpha_a(0.700) < log2/4 < alpha_a(0.701)",
       (lo < alpha_b) and (alpha_b < hi), alpha_a_0700=lo, alpha_a_0701=hi)
# bisection for a tighter enclosure of lambda_x (alpha_a is nondecreasing, Lemma 0)
a_, b_ = arb("0.700"), arb("0.701")
for _ in range(40):
    m = ((a_ + b_) / 2).mid()
    v = alpha_a_arb(m)
    if v < alpha_b:
        a_ = arb(m)
    elif v > alpha_b:
        b_ = arb(m)
    else:
        break
record("lambda_x enclosure (interval comparisons at both ends)",
       (alpha_a_arb(a_.mid()) < alpha_b) and (alpha_b < alpha_a_arb(b_.mid())),
       lo=a_.mid(), hi=b_.mid())

for lam_s in ("0.1", "0.5", "0.7"):
    lam = arb(lam_s)
    c = lam * (1 - lam / 2)
    x = alpha_a_arb(lam_s) / c
    record(f"Lambert root solves I(x) = c x with x < 1 (lambda={lam_s})",
           abs(I(x) - c * x) < arb(10) ** -50 and x < 1, x=x)


# ---------- Part 2: cross-check against the definition (floats, grids) ----------
def A_num(alpha, lam, s_grid, ny):
    """LM eq. (9): inf_s [I(s) + inf_{phi > s} I_{alpha,s,lam}(phi)], with
    I_{alpha,s,lam}(phi) = +inf below phi_typ(s), 0 at phi_typ(s) and
    phi^2/(2s) - alpha above it (LM eq. 64). So the inner infimum is 0 where
    phi_typ(s) > s and s/2 - alpha elsewhere (the code's <= differs only at
    the single point s = s0, which does not change the infimum over s).
    phi_typ(s) (REM free energy) is found by grid maximization of
    (1/lam)(alpha - y^2/(2s)) + y over 0 <= y <= sqrt(2 alpha s); it is
    sqrt(2 alpha s) at lam = inf."""
    ymax = np.sqrt(2 * alpha * s_grid)
    if math.isinf(lam):
        pt = ymax
    else:
        y = ymax[:, None] * np.linspace(0.0, 1.0, ny)[None, :]
        pt = np.max((alpha - y ** 2 / (2 * s_grid[:, None])) / lam + y, axis=1)
    Is = (s_grid - 1 - np.log(s_grid)) / 2
    inner = np.where(s_grid <= pt, 0.0, np.maximum(s_grid / 2 - alpha, 0.0))
    return float(np.min(Is + inner))


def lb_num(lam, s_grid, ny):
    lo_, hi_ = 1e-4, 0.5
    for _ in range(22 if QUICK else 34):
        mid = (lo_ + hi_) / 2
        if A_num(mid, lam, s_grid, ny) - mid > 0:
            lo_ = mid
        else:
            hi_ = mid
    return (lo_ + hi_) / 2


def lb_closed(lam):
    ab = math.log(2) / 4
    if math.isinf(lam) or lam >= 1:
        return min(1 / (2 * math.e), ab)
    return min(float(alpha_a_arb(repr(lam)).mid()), ab)


if QUICK:
    s_grid, ny, tol = np.linspace(0.02, 1.5, 400), 401, 5e-3
    lams = [0.1, 0.5, 0.8, 2.0, math.inf]
else:
    s_grid, ny, tol = np.linspace(0.02, 1.5, 20000), 2001, 2e-4
    lams = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.69, 0.7, 0.705, 0.71, 0.75,
            0.8, 0.9, 0.95, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0, math.inf]
maxdiff = 0.0
for lam in lams:
    n, cf = lb_num(lam, s_grid, ny), lb_closed(lam)
    maxdiff = max(maxdiff, abs(n - cf))
    results["crosscheck"][str(lam)] = {"lb_from_definition": n, "closed_form": cf}
    print(f"lambda={lam}: lb(definition, grid)={n:.6f}  min(alpha_a, log2/4)={cf:.6f}")
cc_pass = maxdiff < tol
ok = ok and cc_pass
results["crosscheck"]["max_abs_diff"] = maxdiff
results["crosscheck"]["tolerance"] = tol
results["crosscheck"]["pass"] = cc_pass
print(f"[{'PASS' if cc_pass else 'FAIL'}] cross-check max |diff| = {maxdiff:.2e} (tol {tol})")

# ---------- Part 3: brute-force support at lambda = inf (full mode only) ----------
if not QUICK:
    rng = np.random.default_rng(SEED)

    def all_retrieved(N, P):
        X = rng.standard_normal((P, N))
        G = X @ X.T
        d = np.diag(G).copy()
        np.fill_diagonal(G, -np.inf)
        return bool(np.all(d > G.max(axis=1)))

    def wilson(k, n, z=1.96):
        p = k / n
        den = 1 + z * z / n
        cen = (p + z * z / (2 * n)) / den
        half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
        return [cen - half, cen + half]

    trials = 200
    for N in (16, 24, 32, 40):
        for alpha in (0.10, 0.17, 0.25):
            if N * alpha > 8.1:  # keep P below about 3000
                continue
            P = int(math.floor(math.exp(alpha * N)))
            k = sum(all_retrieved(N, P) for _ in range(trials))
            results["simulation"][f"N={N},alpha={alpha}"] = {
                "P": P, "trials": trials, "all_retrieved": k, "wilson95": wilson(k, trials)}
            print(f"sim lambda=inf N={N} alpha={alpha} P={P}: {k}/{trials} all retrieved")

results["pass"] = ok
results["mode"] = "quick" if QUICK else "full"
results["seed"] = SEED
results["runtime_s"] = round(time.time() - t_start, 2)
results["versions"] = {"python": platform.python_version(), "numpy": np.__version__,
                       "python-flint": flint.__version__, "arb_precision_bits": ctx.prec}
if not QUICK:
    with open("check-output.json", "w") as f:
        json.dump(results, f, indent=2)
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
