#!/usr/bin/env python
"""Interval-arithmetic (python-flint / Arb) enclosure of the Claim 2.6 Hessian
grad^2 S̄*(1,0) for general kappa, via the general-kappa display in the proof of Lemma B.9
(arXiv:2404.18902 v2; see DEFS.md section 6).

THE PARAMETERS (alpha*, q0, psi0) ARE TAKEN AS INPUT BALLS. The result is a rigorous enclosure
of the Hessian for EVERY (alpha,q,psi) in the box, *provided the integrals are as written*.
That the true (alpha*,q0,psi0) lie in the box is NOT proved here (needs Condition 3.1 certification).

Usage: python flint_hess.py KAPPA ALPHA RA Q RQ PSI RPSI
"""
import sys
from flint import arb, acb, ctx

ctx.prec = 128


def ball(mid, rad):
    return arb(mid, rad)


kappa = arb(sys.argv[1])
al = ball(sys.argv[2], sys.argv[3])
q = ball(sys.argv[4], sys.argv[5])
psi = ball(sys.argv[6], sys.argv[7])
L = 14
sqrt2pi = (2 * arb.pi()).sqrt()


def phi(z):
    return (-z * z / 2).exp() / sqrt2pi


sq, spsi, s1 = q.sqrt(), psi.sqrt(), (1 - q).sqrt()


def gint(f):
    """E f(Z) over [-L,L] (rigorous) -- tail added separately."""
    g = lambda z, an: f(z) * phi(z)
    v = acb.integral(g, -L, L, rel_tol=2.0 ** -60)
    return v.real


def hdot_parts(z):
    t = (spsi * z).tanh()
    om = 1 - t * t
    return t, om, spsi * z


def hhat_parts(z):
    u = (kappa - sq * z) / s1
    E = (2 / arb.pi()).sqrt() * (-u * u / 2).exp() / (u / arb(2).sqrt()).erfc()
    N = E / s1
    Fp = -E * (E - u) / (1 - q)
    Hh = sq * z
    K = u / s1 - N
    return N, Fp, Hh, K


I1 = gint(lambda z: (lambda t, om, h: om * h * h)(*hdot_parts(z)))
I2 = gint(lambda z: (lambda t, om, h: om * h * t)(*hdot_parts(z)))
I3 = gint(lambda z: (lambda t, om, h: om * t * t)(*hdot_parts(z)))
J1 = gint(lambda z: (lambda N, Fp, Hh, K: Fp * (K - Hh / q) ** 2)(*hhat_parts(z)))
J2 = gint(lambda z: (lambda N, Fp, Hh, K: Fp * (K - Hh / q) * N)(*hhat_parts(z)))
J3 = gint(lambda z: (lambda N, Fp, Hh, K: Fp * N * N)(*hhat_parts(z)))
J4 = gint(lambda z: (lambda N, Fp, Hh, K: N * Hh)(*hhat_parts(z)))
J5 = gint(lambda z: (lambda N, Fp, Hh, K: N * N)(*hhat_parts(z)))
J6 = gint(lambda z: (lambda N, Fp, Hh, K: N * K)(*hhat_parts(z)))
# Tail |Z|>L.  Using Lemma 4.21 (0<=E(x)<=|x|+1, 0<E'<1) and q<=0.75, psi<=8, |kappa|<=1, every
# integrand is bounded by C (1+|Z|)^2 with C = 1000 (crude; see REPORT.md).  For |Z|>=1,
# (1+|Z|)^2 <= 4 Z^2 and E[Z^2; |Z|>L] = 2 (L phi(L) + Psi(L)) <= 2 (L+1) phi(L).
assert float(q.lower()) >= 0.5 and float(q.upper()) <= 0.75 and float(psi.upper()) <= 8 and abs(float(kappa)) <= 1
tail = arb(1000) * 4 * 2 * (L + 1) * phi(arb(L))
tball = arb(0, float(tail.upper()) * 1.01)
I1, I2, I3, J1, J2, J3, J4, J5, J6 = [x + tball for x in (I1, I2, I3, J1, J2, J3, J4, J5, J6)]

A = [I2, I3]
B = [I1, I2]
V = [[I1, I2], [I2, I3]]
c = 3 / (1 - q) + 1 / q
Hm = [[None, None], [None, None]]
for i in range(2):
    for j in range(2):
        t1 = -V[i][j]
        t2 = al * (A[i] * A[j] * J1 - (A[i] * B[j] + A[j] * B[i]) / psi * J2 + B[i] * B[j] / psi ** 2 * J3)
        t3 = -al * (-2 / (q * (1 - q)) * A[i] * A[j] * J4 - (A[i] * B[j] + A[j] * B[i]) / (psi * (1 - q)) * J5
                    + c * A[i] * A[j] * J6)
        Hm[i][j] = t1 + t2 + t3
H11, H12, H22 = Hm[0][0], Hm[0][1], Hm[1][1]
det = H11 * H22 - H12 * H12
tr = H11 + H22
disc = (tr * tr - 4 * det)
lam_max = (tr + disc.sqrt()) / 2
lam_min = (tr - disc.sqrt()) / 2
print("tail bound", tail)
for nm, x in [("I1", I1), ("I2", I2), ("I3", I3), ("J1", J1), ("J2", J2), ("J3", J3), ("J4", J4), ("J5", J5), ("J6", J6)]:
    print(nm, x.str(12, radius=True))
print("H11", H11.str(12, radius=True)); print("H12", H12.str(12, radius=True)); print("H22", H22.str(12, radius=True))
print("det", det.str(12, radius=True)); print("lam_max", lam_max.str(12, radius=True)); print("lam_min", lam_min.str(12, radius=True))
print("CERTIFIED negative definite on box:", bool(H11 < 0) and bool(det > 0))
