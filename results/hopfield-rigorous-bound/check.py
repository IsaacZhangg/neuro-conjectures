"""Certify the constant in Theorem 1 of PROOF.md.

The decision step uses only the Python standard library: exact rational
arithmetic (fractions.Fraction) with explicit tail bounds for every series,
so the printed interval for A(tau0) is a proof, not a float estimate.

Optional extras, never part of the pass/fail decision:
- an independent cross-check of the same interval with python-flint (arb);
- a float scan for sup A(tau) and a comparison with Stojnic's closed form;
- an exploratory simulation of zero-temperature dynamics (needs numpy).

Exit 0 when the certified lower bound on A(tau0) is above 0.113, else 1.
Writes check-output.json unless --quick.
"""
import argparse
import json
import math
import sys
from fractions import Fraction as Fr
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TAU0 = Fr(18259, 10000)
THRESHOLD = Fr(113, 1000)
DIGITS = 40  # working precision for square roots, in decimal digits


# ---------- rigorous rational arithmetic ----------

def series_with_tail(term0, ratio, max_terms=400, tol=Fr(1, 10**45)):
    """Sum t_0 + t_1 + ... where t_{n+1} = t_n * ratio(n).

    ratio(n) must be a Fraction whose absolute value is non-increasing in n
    for n past the point where we stop (checked by the callers' formulas:
    every ratio used below is c/(n+1) or c(2n+1)/((n+1)(2n+3)) with c > 0,
    which decreases in n). Returns (s_lo, s_hi) bracketing the full sum,
    using |tail| <= |t_N| / (1 - |r_N|) once |r_N| < 1.
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
    """Interval for exp(-x), x a nonnegative Fraction."""
    return series_with_tail(Fr(1), lambda n: -x / (n + 1))


def atan_inv(m):
    """Interval for arctan(1/m), m a positive integer > 1."""
    x2 = Fr(1, m * m)
    # arctan(x) = sum (-1)^n x^(2n+1)/(2n+1); write t_n = (-1)^n x^(2n+1)/(2n+1)
    return series_with_tail(Fr(1, m), lambda n: -x2 * (2 * n + 1) / (2 * n + 3))


def pi_interval():
    a_lo, a_hi = atan_inv(5)
    b_lo, b_hi = atan_inv(239)
    return 16 * a_lo - 4 * b_hi, 16 * a_hi - 4 * b_lo  # Machin


def sqrt_interval(q):
    """Interval for sqrt(q), q a positive Fraction."""
    scale = 10 ** DIGITS
    n = q.numerator * scale * scale // q.denominator  # floor(q * scale^2)
    lo = Fr(isqrt(n), scale)
    hi = Fr(isqrt(n + 1) + 1, scale)  # sqrt(q*scale^2) < sqrt(n+1) <= isqrt(n+1)+1
    assert lo * lo <= q <= hi * hi
    return lo, hi


def certify(tau):
    """Rigorous intervals for phi(tau), delta(tau)=Q(tau) and A(tau)."""
    pi_lo, pi_hi = pi_interval()
    s_lo, _ = sqrt_interval(2 * pi_lo)
    _, s_hi = sqrt_interval(2 * pi_hi)
    inv_lo, inv_hi = 1 / s_hi, 1 / s_lo  # bounds on 1/sqrt(2 pi)

    e_lo, e_hi = exp_neg(tau * tau / 2)
    phi_lo, phi_hi = inv_lo * e_lo, inv_hi * e_hi

    # Integral_0^tau e^{-u^2/2} du = sum (-1)^n tau^(2n+1) / (2^n n! (2n+1))
    t2 = tau * tau
    i_lo, i_hi = series_with_tail(
        tau, lambda n: -t2 * (2 * n + 1) / (2 * (n + 1) * (2 * n + 3)))
    assert i_lo > 0
    # Q(tau) = 1/2 - (1/sqrt(2 pi)) * integral
    d_lo = Fr(1, 2) - inv_hi * i_hi
    d_hi = Fr(1, 2) - inv_lo * i_lo
    assert 0 < d_lo < d_hi < Fr(1, 2)

    # A = ((delta(1-delta) - phi^2) / phi)^2. delta(1-delta) increases on
    # (0, 1/2); the ratio decreases in phi for positive numerator.
    num_lo = d_lo * (1 - d_lo) - phi_hi * phi_hi
    num_hi = d_hi * (1 - d_hi) - phi_lo * phi_lo
    if num_lo <= 0:
        raise RuntimeError("numerator not certified positive")
    a_lo = (num_lo / phi_hi) ** 2
    a_hi = (num_hi / phi_lo) ** 2
    return {
        "pi": (pi_lo, pi_hi),
        "phi": (phi_lo, phi_hi),
        "delta": (d_lo, d_hi),
        "A": (a_lo, a_hi),
    }


def dec(x, places=12):
    """Decimal string of a Fraction, rounded toward minus infinity."""
    return _fmt(math.floor(x * 10**places), places)


def dec_up(x, places=12):
    """Decimal string of a Fraction, rounded toward plus infinity."""
    return _fmt(math.ceil(x * 10**places), places)


def _fmt(n, places):
    sign = "-" if n < 0 else ""
    n = abs(n)
    scale = 10 ** places
    return f"{sign}{n // scale}.{n % scale:0{places}d}"


# ---------- optional cross-checks ----------

def arb_crosscheck(tau):
    try:
        import flint
    except ImportError:
        return {"ran": False, "reason": "python-flint not installed"}
    flint.ctx.prec = 200
    t = flint.arb(tau.numerator) / tau.denominator
    phi = (-(t * t) / 2).exp() / (2 * flint.arb.pi()).sqrt()
    delta = (t / flint.arb(2).sqrt()).erfc() / 2
    a = ((delta * (1 - delta) - phi * phi) / phi) ** 2
    lo = float(a.lower())
    hi = float(a.upper())
    return {"ran": True, "A_lower": lo, "A_upper": hi, "agrees": True}


def float_scan():
    def a(t):
        phi = math.exp(-t * t / 2) / math.sqrt(2 * math.pi)
        d = 0.5 * math.erfc(t / math.sqrt(2))
        num = d * (1 - d) - phi * phi
        return (num / phi) ** 2 if num > 0 else 0.0

    lo, hi = 1.5, 2.2
    for _ in range(200):  # golden-section search, A is unimodal here
        m1 = lo + (hi - lo) * 0.381966
        m2 = lo + (hi - lo) * 0.618034
        if a(m1) < a(m2):
            lo = m1
        else:
            hi = m2
    t = (lo + hi) / 2

    # Stojnic (2024) reports alpha_c(NLT,1) = erf(x)^2/(2x^2) - 1 + erf(x)^2
    # at x = 1.29108. Evaluate his expression at x = tau*/sqrt(2) and compare.
    x = t / math.sqrt(2)
    e = math.erf(x)
    stojnic = e * e / (2 * x * x) - 1 + e * e
    return {
        "tau_star": round(t, 10),
        "delta_star": round(0.5 * math.erfc(t / math.sqrt(2)), 10),
        "sup_A_float": round(a(t), 10),
        "stojnic_x": round(x, 8),
        "stojnic_expression_at_x": round(stojnic, 10),
        "note": "floating point, not certified; sup A is not used by the theorem",
    }


def simulate(quick, seed):
    """Exploratory: zero-temperature asynchronous dynamics started at the
    target pattern, for model G (Gaussian interference) and model R (all
    patterns +-1). Reports the final overlap. Not part of the decision."""
    try:
        import numpy as np
    except ImportError:
        return {"ran": False, "reason": "numpy not installed"}
    rng = np.random.default_rng(seed)
    n = 300 if quick else 1000
    trials = 3 if quick else 10
    rows = []
    for model in ("G", "R"):
        for alpha in (0.10, 0.12, 0.129):
            p = int(round(alpha * n))
            overlaps = []
            for _ in range(trials):
                xi = np.ones(n)
                if model == "G":
                    other = rng.standard_normal((p - 1, n))
                else:
                    other = rng.choice([-1.0, 1.0], size=(p - 1, n))
                pats = np.vstack([xi, other])
                J = pats.T @ pats / n
                np.fill_diagonal(J, 0.0)
                s = xi.copy()
                h = J @ s
                for _sweep in range(100):
                    changed = 0
                    for i in rng.permutation(n):
                        if s[i] * h[i] < 0:
                            s[i] = -s[i]
                            h += 2 * s[i] * J[:, i]
                            changed += 1
                    if changed == 0:
                        break
                overlaps.append(float(s @ xi / n))
            rows.append({"model": model, "alpha": alpha, "N": n,
                         "trials": trials,
                         "mean_final_overlap": round(sum(overlaps) / trials, 4),
                         "min_final_overlap": round(min(overlaps), 4)})
    return {"ran": True, "rows": rows,
            "note": "exploratory evidence only; finite N, not a proof"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--seed", type=int, default=20261125)
    args = parser.parse_args()

    c = certify(TAU0)
    a_lo, a_hi = c["A"]
    passed = a_lo > THRESHOLD
    print(f"tau0 = {TAU0} = {float(TAU0)}")
    print(f"delta(tau0) in [{dec(c['delta'][0])}, {dec_up(c['delta'][1])}]")
    print(f"phi(tau0)   in [{dec(c['phi'][0])}, {dec_up(c['phi'][1])}]")
    print(f"A(tau0)     in [{dec(a_lo)}, {dec_up(a_hi)}]  (certified, exact rationals)")
    print(f"certified A(tau0) > {float(THRESHOLD)}: {passed}")
    print(f"certified A(tau0) > 0.12948: {a_lo > Fr(12948, 100000)}")

    arb = arb_crosscheck(TAU0)
    if arb["ran"]:
        arb["agrees"] = arb["A_lower"] <= float(a_hi) and float(a_lo) <= arb["A_upper"]
        print(f"arb cross-check: [{arb['A_lower']:.12f}, {arb['A_upper']:.12f}] agrees={arb['agrees']}")
    scan = float_scan()
    print(f"float sup A = {scan['sup_A_float']} at tau = {scan['tau_star']}; "
          f"Stojnic's expression at x = tau*/sqrt(2) = {scan['stojnic_expression_at_x']}")
    # Theorem 3 (reference-state refinement), see PROOF-refstate.md.
    import check_refstate as cr
    out3, checks3 = cr.certify(cr.TAU0, cr.THETA0)
    passed3 = all(bool(v) for v in checks3.values())
    a2_lo, a2_hi = out3["A2"]
    print(f"Theorem 3: A2({cr.TAU0}, {cr.THETA0}) in [{dec(a2_lo)}, {dec_up(a2_hi)}]"
          f"  (certified); side conditions hold: {passed3}")
    arb3 = cr.arb_crosscheck(cr.TAU0, cr.THETA0)
    if arb3["ran"]:
        arb3 = {"ran": True, "A2_lower": float(arb3["A2_lower"]),
                "A2_upper": float(arb3["A2_upper"])}
        arb3["agrees"] = (arb3["A2_lower"] <= float(a2_hi) + 1e-15
                          and float(a2_lo) <= arb3["A2_upper"] + 1e-15)
        print(f"Theorem 3 arb cross-check agrees: {arb3['agrees']}")
    sim = simulate(args.quick, args.seed)
    if sim["ran"]:
        for r in sim["rows"]:
            print(f"  sim model {r['model']} alpha={r['alpha']}: mean overlap "
                  f"{r['mean_final_overlap']}, min {r['min_final_overlap']}")

    if not args.quick:
        out = {
            "result": "hopfield-rigorous-bound",
            "mode": "full",
            "seed": args.seed,
            "decision": {
                "rule": "certified lower bound on A(tau0) exceeds 0.113",
                "tau0": str(TAU0),
                "A_tau0_lower": dec(a_lo, 15),
                "A_tau0_upper": dec_up(a_hi, 15),
                "delta_tau0": [dec(c["delta"][0], 15), dec_up(c["delta"][1], 15)],
                "phi_tau0": [dec(c["phi"][0], 15), dec_up(c["phi"][1], 15)],
                "passed": passed,
                "method": "exact rational arithmetic with explicit series tail bounds",
            },
            "arb_crosscheck": arb,
            "theorem3": {
                "rule": "certified lower bound on A2(tau0, theta0) with side conditions (reported; added after the original decision rule)",
                "tau0": str(cr.TAU0),
                "theta0": str(cr.THETA0),
                "A2_lower": dec(a2_lo, 15),
                "A2_upper": dec_up(a2_hi, 15),
                "checks": {k: bool(v) for k, v in checks3.items()},
                "passed": passed3,
                "arb_crosscheck": arb3,
            },
            "float_scan": scan,
            "simulation": sim,
            "scope": ("Theorem 1 (Gaussian interference) is proved; the +-1 "
                      "Hopfield statement (Theorem 2) is conditional on "
                      "hypothesis (U) and is not proved."),
        }
        (HERE / "check-output.json").write_text(json.dumps(out, indent=2) + "\n")
    ok = passed and (not arb["ran"] or arb["agrees"])
    ok = ok and passed3 and (not arb3["ran"] or arb3["agrees"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
