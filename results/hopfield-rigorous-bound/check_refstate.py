"""Certify the constant in Theorem 3 of PROOF-refstate.md.

Decision step: Python standard library only. Exact rational arithmetic
(fractions.Fraction) with explicit tail bounds for every series, the same
method as check.py of the original hopfield-rigorous-bound result, extended to
phi and Q at the second point theta.

Certified (as interval statements with rational endpoints):
  * 0 < tau0 < theta0, so Q(theta0) < Q(tau0) is also checked from intervals;
  * D   = phi(tau0) - phi(theta0)(1 - 2Q(theta0))                > 0;
  * Num = Q(tau0)(1-Q(tau0)) - phi(tau0)^2
          + phi(theta0)^2 - Q(theta0)(1 - Q(theta0))             > 0;
  * A2(tau0, theta0) = (Num/D)^2 >= THRESHOLD.

Optional, never part of the decision: an arb (python-flint) cross-check and a
float optimisation of A2 over (tau, theta).

Exit 0 if every certified inequality holds (and the arb cross-check, when it
runs, agrees), else 1.
"""
import math
import sys
from fractions import Fraction as Fr
from math import isqrt

TAU0 = Fr(9711, 5000)      # 1.9422
THETA0 = Fr(1582, 625)     # 2.5312
THRESHOLD = Fr(134716, 1000000)
DIGITS = 40                # decimal digits for rational square-root brackets


# ---------- rigorous rational arithmetic (as in check.py) ----------

def series_with_tail(term0, ratio, max_terms=600, tol=Fr(1, 10**45)):
    """Bracket t_0 + t_1 + ... where t_{n+1} = t_n * ratio(n).

    Every ratio used below is -x/(n+1) or -c(2n+1)/((n+1)(2n+3)) with
    x, c > 0, whose absolute value is non-increasing in n. Once |r_n| < 1 the
    remainder sum_{j>=n} t_j has absolute value at most |t_n|/(1 - |r_n|),
    because |t_{j+1}| = |t_j||r_j| <= |t_j||r_n| for j >= n.
    """
    s = Fr(0)
    t = term0
    for n in range(max_terms):
        r = ratio(n)
        if abs(r) < 1 and abs(t) / (1 - abs(r)) < tol:
            bound = abs(t) / (1 - abs(r))
            return s - bound, s + bound
        s += t
        t = t * r
    raise RuntimeError("series did not converge to tolerance")


def exp_neg(x):
    """Interval for exp(-x), x >= 0 rational."""
    return series_with_tail(Fr(1), lambda n: -x / (n + 1))


def atan_inv(m):
    """Interval for arctan(1/m), m > 1 an integer."""
    x2 = Fr(1, m * m)
    return series_with_tail(Fr(1, m), lambda n: -x2 * (2 * n + 1) / (2 * n + 3))


def pi_interval():
    a_lo, a_hi = atan_inv(5)
    b_lo, b_hi = atan_inv(239)
    return 16 * a_lo - 4 * b_hi, 16 * a_hi - 4 * b_lo  # Machin's formula


def sqrt_interval(q):
    """Interval for sqrt(q), q > 0 rational."""
    scale = 10 ** DIGITS
    n = q.numerator * scale * scale // q.denominator
    lo = Fr(isqrt(n), scale)
    hi = Fr(isqrt(n + 1) + 1, scale)
    assert lo * lo <= q <= hi * hi
    return lo, hi


def inv_sqrt_2pi():
    pi_lo, pi_hi = pi_interval()
    s_lo, _ = sqrt_interval(2 * pi_lo)
    _, s_hi = sqrt_interval(2 * pi_hi)
    return 1 / s_hi, 1 / s_lo


def phi_Q(t, inv):
    """Rigorous intervals for phi(t) and Q(t) = 1 - Phi(t), t > 0 rational."""
    inv_lo, inv_hi = inv
    e_lo, e_hi = exp_neg(t * t / 2)
    assert e_lo > 0
    phi = (inv_lo * e_lo, inv_hi * e_hi)
    # int_0^t e^{-u^2/2} du = sum_n (-1)^n t^(2n+1) / (2^n n! (2n+1))
    t2 = t * t
    i_lo, i_hi = series_with_tail(
        t, lambda n: -t2 * (2 * n + 1) / (2 * (n + 1) * (2 * n + 3)))
    assert i_lo > 0
    Q = (Fr(1, 2) - inv_hi * i_hi, Fr(1, 2) - inv_lo * i_lo)
    assert 0 < Q[0] <= Q[1] < Fr(1, 2)
    return phi, Q


def certify(tau, theta):
    inv = inv_sqrt_2pi()
    (pt_lo, pt_hi), (d_lo, d_hi) = phi_Q(tau, inv)
    (pq_lo, pq_hi), (q_lo, q_hi) = phi_Q(theta, inv)
    out = {"phi_tau": (pt_lo, pt_hi), "Q_tau": (d_lo, d_hi),
           "phi_theta": (pq_lo, pq_hi), "Q_theta": (q_lo, q_hi)}
    checks = {}
    checks["0 < tau0 < theta0"] = 0 < tau < theta
    checks["Q(theta0) < Q(tau0)"] = q_hi < d_lo
    checks["1 - 2Q(theta0) > 0"] = 1 - 2 * q_hi > 0
    # x -> x(1-x) is increasing on (0, 1/2); both Q's lie in (0, 1/2).
    num_lo = d_lo * (1 - d_lo) - pt_hi ** 2 + pq_lo ** 2 - q_hi * (1 - q_hi)
    num_hi = d_hi * (1 - d_hi) - pt_lo ** 2 + pq_hi ** 2 - q_lo * (1 - q_lo)
    # D = phi_tau - phi_theta (1 - 2 Q_theta); phi_theta >= 0, 1-2Q_theta > 0.
    den_lo = pt_lo - pq_hi * (1 - 2 * q_lo)
    den_hi = pt_hi - pq_lo * (1 - 2 * q_hi)
    checks["Num > 0"] = num_lo > 0
    checks["D > 0"] = den_lo > 0
    out["Num"] = (num_lo, num_hi)
    out["D"] = (den_lo, den_hi)
    if num_lo > 0 and den_lo > 0:
        out["A2"] = ((num_lo / den_hi) ** 2, (num_hi / den_lo) ** 2)
        checks[f"A2(tau0, theta0) >= {THRESHOLD}"] = out["A2"][0] >= THRESHOLD
    else:
        checks[f"A2(tau0, theta0) >= {THRESHOLD}"] = False
    return out, checks


def dec(x, places=12):
    return _fmt(math.floor(x * 10**places), places)


def dec_up(x, places=12):
    return _fmt(math.ceil(x * 10**places), places)


def _fmt(n, places):
    sign = "-" if n < 0 else ""
    n = abs(n)
    scale = 10 ** places
    return f"{sign}{n // scale}.{n % scale:0{places}d}"


# ---------- optional extras (not part of the decision) ----------

def arb_crosscheck(tau, theta):
    try:
        import flint
    except ImportError:
        return {"ran": False, "reason": "python-flint not installed"}
    flint.ctx.prec = 256

    def pq(x):
        t = flint.arb(x.numerator) / x.denominator
        return ((-(t * t) / 2).exp() / (2 * flint.arb.pi()).sqrt(),
                (t / flint.arb(2).sqrt()).erfc() / 2)

    ft, d = pq(tau)
    fq, q = pq(theta)
    num = d * (1 - d) - ft * ft + fq * fq - q * (1 - q)
    den = ft - fq * (1 - 2 * q)
    a2 = (num / den) ** 2
    return {"ran": True, "A2_lower": a2.lower(), "A2_upper": a2.upper(),
            "Q_theta_upper": q.upper(), "Q_tau_lower": d.lower()}


def float_optimum():
    def phi(t):
        return math.exp(-t * t / 2) / math.sqrt(2 * math.pi)

    def Q(t):
        return 0.5 * math.erfc(t / math.sqrt(2))

    def a2(t, th):
        if not 0 < t < th:
            return 0.0
        num = Q(t) * (1 - Q(t)) - phi(t) ** 2 + phi(th) ** 2 - Q(th) * (1 - Q(th))
        den = phi(t) - phi(th) * (1 - 2 * Q(th))
        return (num / den) ** 2 if num > 0 and den > 0 else 0.0

    best = (0.0, 0.0, 0.0)
    for i in range(121):
        for j in range(121):
            t = 1.7 + 0.004 * i
            th = 2.3 + 0.005 * j
            v = a2(t, th)
            if v > best[0]:
                best = (v, t, th)
    v, t, th = best
    h = 0.004
    for _ in range(60):  # crude coordinate refinement
        for dt, dth in ((h, 0), (-h, 0), (0, h), (0, -h)):
            w = a2(t + dt, th + dth)
            if w > v:
                v, t, th = w, t + dt, th + dth
        h *= 0.7
    return v, t, th


def main() -> int:
    out, checks = certify(TAU0, THETA0)
    print(f"tau0   = {TAU0} = {float(TAU0)}")
    print(f"theta0 = {THETA0} = {float(THETA0)}")
    for key in ("phi_tau", "Q_tau", "phi_theta", "Q_theta", "Num", "D"):
        lo, hi = out[key]
        print(f"{key:10s} in [{dec(lo, 15)}, {dec_up(hi, 15)}]")
    if "A2" in out:
        lo, hi = out["A2"]
        print(f"A2(tau0,theta0) in [{dec(lo, 12)}, {dec_up(hi, 12)}]  "
              "(certified, exact rationals)")
    ok = True
    for name, val in checks.items():
        print(f"certified {name}: {val}")
        ok = ok and bool(val)

    arb = arb_crosscheck(TAU0, THETA0)
    if arb["ran"]:
        lo, hi = out["A2"]
        agrees = (float(arb["A2_lower"]) <= float(hi) + 1e-15
                  and float(lo) <= float(arb["A2_upper"]) + 1e-15)
        print(f"arb cross-check A2 in [{float(arb['A2_lower']):.15f}, {float(arb['A2_upper']):.15f}] "
              f"agrees={agrees}")
        ok = ok and agrees
    else:
        print(f"arb cross-check skipped: {arb['reason']}")
    v, t, th = float_optimum()
    print(f"float (uncertified) sup A2 ~ {v:.7f} at tau ~ {t:.4f}, theta ~ {th:.4f}")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
