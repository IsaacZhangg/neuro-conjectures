"""Certify every numeric claim used by PROOF-pm1.md (Theorem 2, +-1 patterns).

Decision step: Python standard library only, exact rational arithmetic
(fractions.Fraction). Every transcendental quantity (pi, exp, the normal
integral, ln) is enclosed in an interval with rational endpoints using a
series plus an explicit tail bound; square roots use integer-sqrt brackets.
The helpers series_with_tail, exp_neg, atan_inv, pi_interval and
sqrt_interval are copied from check.py of hopfield-rigorous-bound (same
method). ln_interval is new; its tail bound is stated in its docstring.

Certified claims (labels match PROOF-pm1.md, section 8):
  C1  0 < delta0 = Q(tau0) < 1/2 and L0 = ln(1/delta0) > 1/2 (so t_H >= 0),
      and tau0 <= sqrt(2 L0) (so every t used is in [0, sqrt(2 L0)]).
  C2  D(tau0) = delta0*sqrt(2 L0) - phi(tau0) lies in a certified interval
      with D > 0.
  C3  condition (a) at alpha = 3/100, rho0 = 17/20:  2D/(1-rho0) <= sqrt(alpha)
      (checked as (2D/(1-rho0))^2 <= alpha), hence for all alpha >= 3/100.
  C4  condition (b) at alpha = 3/100: 2*phi(0)(1-rho0)(1+rho0)^2/(6 rho0^3)
      <= sqrt(alpha) (checked squared), hence for all alpha >= 3/100.
  C5  A(tau0) > 12948/100000.
  C6  A_H(delta0) = (1 - delta0 - 2 delta0 L0)^2 / (2 L0) > 3/100
      (Hoeffding route at the SAME radius delta0; this is what the proof uses).
  C7  A_H at the rational tau_H = 22217/10000: A_H(Q(tau_H)) > 0.0879 (> 0.03).
      Informational: the best Hoeffding-only constant; not needed by the proof.
  C8  analytic constants used in the proof:
      2 e^{-3/2} <= 1                      (sup|phi''| = phi(0));
      2 sqrt(2 ln 9) <= 21/5               (operator-norm lemma);
      (1 + 2 sqrt(2/pi))/3 <= 1            (Lindeberg third-moment constant);
      8 sqrt(pi) + 2/sqrt(ln 4) <= 16      (Step 2 constant);
      8 sqrt(pi) + 4/sqrt(ln 5) <= 18      (Step 5 constant);
      1 - 2 delta0 >= 0.932                (overlap).
  C9  the barrier inequality itself at the two ends used, as a consistency
      check of the algebra: (sqrt(a) + 2 phi(tau0))^2 < a + 4 delta0(1-delta0)
      at a = 12948/100000, and (sqrt(a) + 2 phi_H)^2 < a + 4 delta0(1-delta0)
      at a = 3/100 (phi_H = delta0 sqrt(2 L0)).

Optional, never part of the decision: python-flint (arb) cross-check of the
same quantities at 256 bits.

Exit code 0 iff every certified check passes (and the arb cross-check, when
it runs, agrees with the rational intervals).
"""
import sys
from fractions import Fraction as Fr
from math import isqrt, floor, ceil

TAU0 = Fr(18259, 10000)
RHO0 = Fr(17, 20)
ALPHA_LO = Fr(3, 100)          # lower end of Case A / upper end of Case B
ALPHA_HI = Fr(12948, 100000)   # stated upper threshold
TAU_H = Fr(22217, 10000)       # informational Hoeffding optimum
DIGITS = 45


# ---------- rigorous rational arithmetic (copied from check.py) ----------

def series_with_tail(term0, ratio, max_terms=800, tol=Fr(1, 10**50)):
    """Bracket t_0 + t_1 + ... where t_{n+1} = t_n * ratio(n).

    Every ratio passed below is -x/(n+1) or -c(2n+1)/((n+1)(2n+3)) with
    x, c > 0, whose absolute value is non-increasing in n. Once |r_n| < 1 the
    remainder sum_{j>=n} t_j has absolute value at most |t_n|/(1 - |r_n|),
    because |t_{j+1}| <= |t_j||r_n| for j >= n.
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
    assert x >= 0
    return series_with_tail(Fr(1), lambda n: -x / (n + 1))


def atan_inv(m):
    x2 = Fr(1, m * m)
    return series_with_tail(Fr(1, m), lambda n: -x2 * (2 * n + 1) / (2 * n + 3))


def pi_interval():
    a_lo, a_hi = atan_inv(5)
    b_lo, b_hi = atan_inv(239)
    return 16 * a_lo - 4 * b_hi, 16 * a_hi - 4 * b_lo  # Machin


def sqrt_interval(q):
    """Interval [lo, hi] with lo^2 <= q <= hi^2, q > 0 rational."""
    assert q > 0
    scale = 10 ** DIGITS
    n = q.numerator * scale * scale // q.denominator
    lo = Fr(isqrt(n), scale)
    hi = Fr(isqrt(n + 1) + 1, scale)
    assert lo * lo <= q <= hi * hi
    return lo, hi


def sqrt_lo(q):
    return sqrt_interval(q)[0]


def sqrt_hi(q):
    return sqrt_interval(q)[1]


# ---------- new: natural logarithm ----------

def _ln_1_2(y, tol=Fr(1, 10**50)):
    """Interval for ln(y), 1 <= y <= 2 rational.

    ln y = 2 atanh(z) = sum_{n>=0} 2 z^(2n+1)/(2n+1), z = (y-1)/(y+1) in
    [0, 1/3]. All terms are >= 0, so a partial sum is a lower bound; the
    remainder after n terms is at most 2 z^(2n+1)/((2n+1)(1 - z^2)) since
    the term ratio z^2 (2j+1)/(2j+3) is below z^2.
    """
    assert 1 <= y <= 2
    z = (y - 1) / (y + 1)
    z2 = z * z
    s = Fr(0)
    p = z  # z^(2n+1)
    for n in range(2000):
        term = 2 * p / (2 * n + 1)
        tail = term / (1 - z2)
        if tail < tol:
            return s, s + tail
        s += term
        p *= z2
    raise RuntimeError("ln series did not converge")


def ln_interval(q):
    """Interval for ln(q), q > 0 rational."""
    assert q > 0
    if q < 1:
        lo, hi = ln_interval(1 / q)
        return -hi, -lo
    m = 0
    while q >= 2 ** (m + 1):
        m += 1
    y = q / (2 ** m)
    l2_lo, l2_hi = _ln_1_2(Fr(2))
    y_lo, y_hi = _ln_1_2(y)
    return m * l2_lo + y_lo, m * l2_hi + y_hi


# ---------- normal density and tail ----------

def inv_sqrt_2pi():
    pi_lo, pi_hi = pi_interval()
    return 1 / sqrt_hi(2 * pi_hi), 1 / sqrt_lo(2 * pi_lo)


def phi_Q(t):
    """Intervals for phi(t) and Q(t) = P(g > t), t > 0 rational."""
    inv_lo, inv_hi = inv_sqrt_2pi()
    e_lo, e_hi = exp_neg(t * t / 2)
    phi = (inv_lo * e_lo, inv_hi * e_hi)
    t2 = t * t
    i_lo, i_hi = series_with_tail(
        t, lambda n: -t2 * (2 * n + 1) / (2 * (n + 1) * (2 * n + 3)))
    assert i_lo > 0
    Q = (Fr(1, 2) - inv_hi * i_hi, Fr(1, 2) - inv_lo * i_lo)
    assert 0 < Q[0] < Q[1] < Fr(1, 2)
    return phi, Q


def dec(x, places=12):
    return _fmt(floor(x * 10**places), places)


def dec_up(x, places=12):
    return _fmt(ceil(x * 10**places), places)


def _fmt(n, places):
    sign = "-" if n < 0 else ""
    n = abs(n)
    s = 10 ** places
    return f"{sign}{n // s}.{n % s:0{places}d}"


def iv(p, places=12):
    return f"[{dec(p[0], places)}, {dec_up(p[1], places)}]"


# ---------- the certified computation ----------

def hoeffding_AH(d, L):
    """Rigorous lower/upper bounds for A_H = (1 - d - 2 d L)^2/(2 L)."""
    num_lo = 1 - d[1] - 2 * d[1] * L[1]
    num_hi = 1 - d[0] - 2 * d[0] * L[0]
    if num_lo <= 0:
        raise RuntimeError("A_H numerator not certified positive")
    return num_lo * num_lo / (2 * L[1]), num_hi * num_hi / (2 * L[0])


def certify():
    out, checks = {}, {}
    (ph_lo, ph_hi), (d_lo, d_hi) = phi_Q(TAU0)
    out["phi(tau0)"] = (ph_lo, ph_hi)
    out["delta0"] = (d_lo, d_hi)
    phi0 = inv_sqrt_2pi()
    out["phi(0)"] = phi0

    # L0 = ln(1/delta0); ln is increasing, 1/delta decreasing in delta
    L_lo = ln_interval(1 / d_hi)[0]
    L_hi = ln_interval(1 / d_lo)[1]
    out["L0"] = (L_lo, L_hi)
    checks["C1 delta0 in (0,1/2)"] = 0 < d_lo and d_hi < Fr(1, 2)
    checks["C1 L0 > 1/2"] = L_lo > Fr(1, 2)
    checks["C1 tau0^2 <= 2 L0"] = TAU0 * TAU0 <= 2 * L_lo

    # phi_H = delta0 * sqrt(2 L0), D = phi_H - phi(tau0)
    pH = (d_lo * sqrt_lo(2 * L_lo), d_hi * sqrt_hi(2 * L_hi))
    out["phi_H = delta0 sqrt(2 L0)"] = pH
    D = (pH[0] - ph_hi, pH[1] - ph_lo)
    out["D(tau0)"] = D
    checks["C2 D(tau0) > 0"] = D[0] > 0

    # (a): D <= (1-rho0) sqrt(alpha)/2  <=>  (2D/(1-rho0))^2 <= alpha (D>0)
    a_thr = (2 * D[1] / (1 - RHO0)) ** 2
    out["(a) threshold alpha_a <="] = a_thr
    checks["C3 (a) at alpha=3/100"] = a_thr <= ALPHA_LO
    checks["C3 (a) at alpha=12948/100000"] = a_thr <= ALPHA_HI

    # (b): phi(0)(1-rho0)(1+rho0)^2/(6 rho0^3) <= sqrt(alpha)/2
    bL = phi0[1] * (1 - RHO0) * (1 + RHO0) ** 2 / (6 * RHO0 ** 3)
    out["(b) left side <="] = bL
    b_thr = (2 * bL) ** 2
    out["(b) threshold alpha_b <="] = b_thr
    checks["C4 (b) at alpha=3/100"] = b_thr <= ALPHA_LO
    checks["C4 (b) at alpha=12948/100000"] = b_thr <= ALPHA_HI

    # A(tau0) = ((d(1-d) - phi^2)/phi)^2
    num_lo = d_lo * (1 - d_lo) - ph_hi * ph_hi
    num_hi = d_hi * (1 - d_hi) - ph_lo * ph_lo
    assert num_lo > 0
    A = ((num_lo / ph_hi) ** 2, (num_hi / ph_lo) ** 2)
    out["A(tau0)"] = A
    checks["C5 A(tau0) > 0.12948"] = A[0] > ALPHA_HI

    # Hoeffding route at delta0
    AH0 = hoeffding_AH((d_lo, d_hi), (L_lo, L_hi))
    out["A_H(delta0)"] = AH0
    checks["C6 A_H(delta0) > 0.03"] = AH0[0] > ALPHA_LO

    # Hoeffding route at tau_H (informational)
    _, (dh_lo, dh_hi) = phi_Q(TAU_H)
    LH = (ln_interval(1 / dh_hi)[0], ln_interval(1 / dh_lo)[1])
    AHH = hoeffding_AH((dh_lo, dh_hi), LH)
    out["delta_H = Q(tau_H)"] = (dh_lo, dh_hi)
    out["A_H(Q(tau_H))"] = AHH
    checks["C7 A_H(Q(tau_H)) > 0.0879"] = AHH[0] > Fr(879, 10000)
    checks["C7 A_H(Q(tau_H)) > 0.03"] = AHH[0] > ALPHA_LO

    # analytic constants
    e15 = exp_neg(Fr(3, 2))
    checks["C8 2 e^{-3/2} <= 1"] = 2 * e15[1] <= 1
    ln9 = ln_interval(Fr(9))
    checks["C8 2 sqrt(2 ln 9) <= 21/5"] = 4 * 2 * ln9[1] <= Fr(21, 5) ** 2
    pi_lo, pi_hi = pi_interval()
    checks["C8 (1 + 2 sqrt(2/pi))/3 <= 1"] = (1 + 2 * sqrt_hi(2 / pi_lo)) <= 3
    ln4 = ln_interval(Fr(4))
    ln5 = ln_interval(Fr(5))
    s_pi = sqrt_hi(pi_hi)
    checks["C8 8 sqrt(pi) + 2/sqrt(ln 4) <= 16"] = 8 * s_pi + 2 / sqrt_lo(ln4[0]) <= 16
    checks["C8 8 sqrt(pi) + 4/sqrt(ln 5) <= 18"] = 8 * s_pi + 4 / sqrt_lo(ln5[0]) <= 18
    checks["C8 1 - 2 delta0 >= 0.932"] = 1 - 2 * d_hi >= Fr(932, 1000)

    # barrier inequality at the ends (consistency of the algebra)
    rhs_lo_hi = lambda a: a + 4 * d_lo * (1 - d_lo)   # d(1-d) increasing on (0,1/2)
    lhsA = (sqrt_hi(ALPHA_HI) + 2 * ph_hi) ** 2
    checks["C9 barrier Case A at 0.12948"] = lhsA < rhs_lo_hi(ALPHA_HI)
    lhsB = (sqrt_hi(ALPHA_LO) + 2 * pH[1]) ** 2
    checks["C9 barrier Case B at 0.03"] = lhsB < rhs_lo_hi(ALPHA_LO)
    out["C9 margins (rhs_lo - lhs_hi)"] = (rhs_lo_hi(ALPHA_HI) - lhsA,
                                           rhs_lo_hi(ALPHA_LO) - lhsB)
    return out, checks


def arb_crosscheck(out):
    try:
        import flint
    except ImportError:
        return {"ran": False, "reason": "python-flint not installed"}
    flint.ctx.prec = 256
    A = flint.arb

    def fr(x):
        return A(x.numerator) / x.denominator

    def pq(t):
        t = fr(t)
        return (-(t * t) / 2).exp() / (2 * A.pi()).sqrt(), (t / A(2).sqrt()).erfc() / 2

    ph, d = pq(TAU0)
    L = (1 / d).log()
    pH = d * (2 * L).sqrt()
    D = pH - ph
    Aval = ((d * (1 - d) - ph * ph) / ph) ** 2
    AH0 = (1 - d - 2 * d * L) ** 2 / (2 * L)
    _, dh = pq(TAU_H)
    LH = (1 / dh).log()
    AHH = (1 - dh - 2 * dh * LH) ** 2 / (2 * LH)
    vals = {"D(tau0)": D, "A(tau0)": Aval, "A_H(delta0)": AH0,
            "A_H(Q(tau_H))": AHH, "L0": L, "delta0": d}
    agree = True
    rep = {}
    for k, v in vals.items():
        lo, hi = out[k]
        # arb ball must intersect the rational interval
        a_lo, a_hi = v.lower(), v.upper()
        ok = (a_lo <= fr(hi)) and (fr(lo) <= a_hi)
        rep[k] = (v.str(15), bool(ok))
        agree = agree and bool(ok)
    return {"ran": True, "values": rep, "agrees": agree}


def main() -> int:
    out, checks = certify()
    print("check_pm1.py: exact-rational certification for PROOF-pm1.md")
    print(f"tau0 = {TAU0} = {float(TAU0)}, rho0 = {RHO0}, tau_H = {TAU_H}")
    for k, v in out.items():
        if isinstance(v, tuple):
            print(f"  {k:32s} in {iv(v)}")
        else:
            print(f"  {k:32s} {dec_up(v)}")
    print("certified checks:")
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    passed = all(checks.values())
    arb = arb_crosscheck(out)
    if arb["ran"]:
        print("arb cross-check (256 bits):")
        for k, (s, ok) in arb["values"].items():
            print(f"  {k:16s} {s}  agrees={ok}")
    else:
        print("arb cross-check skipped:", arb["reason"])
    ok = passed and (not arb["ran"] or arb["agrees"])
    print("ALL PASS" if ok else "SOME CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
