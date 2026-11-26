"""Certified enclosure of the replica-symmetric fixed point (alpha_*, q_*, psi_*) at margin kappa.

Agent-written (Claude), not independently reviewed.  Reuses the KM class of certify_1d.py
(python-flint / Arb ball arithmetic).  Certifies, for the box
    alpha in [a_lb, a_ub],  q in [q_lb, q_ub]:
  (B1) sup dP(R(q,alpha))/dq < 1 on the box                      (Ding-Sun (5) / AT-type contraction)
  (B2) T(a_lb,q_lb)-q_lb > 0 > T(a_ub,q_ub)-q_ub, T = P o R        (fixed point inside, unique)
  (B3) G_*(a_lb) > 0 > G_*(a_ub) and dG_*/dalpha < 0               (alpha_* in (a_lb,a_ub))
and returns enclosures Q0 of {q_*(alpha)} and Psi0 of {psi_*(alpha)} over alpha in [a_lb,a_ub].
"""
import sys, json, time
sys.path.insert(0, '/home/claude/work')
from certify_1d import KM, box, fixed_point_mid
from arbgauss import ub
from flint import arb, ctx

def fp_box(kappa, alpha_f, q_f, dalpha=1e-9, dq=1e-7, prec=128, verbose=True):
    ctx.prec = prec
    km = KM(kappa)
    t0 = time.time()
    a_lb = arb((arb(alpha_f) - arb(dalpha)).mid()); a_ub = arb((arb(alpha_f) + arb(dalpha)).mid())
    A = box(a_lb, a_ub)
    q_lb = arb(arb(q_f - dq).mid()); q_ub = arb(arb(q_f + dq).mid())
    Q = box(q_lb, q_ub)
    dTmax = km.dT(A, Q)
    c1 = ub(dTmax) < 1
    s_lb = km.T(a_lb, q_lb) - q_lb
    s_ub = km.T(a_ub, q_ub) - q_ub
    c2 = bool((s_lb > 0) and (s_ub < 0))
    def tight(alpha):
        qm = fixed_point_mid(km, alpha, arb(q_f))
        e = arb('1e-28')
        lo = arb((qm - e).mid()); hi = arb((qm + e).mid())
        ok = (km.T(alpha, lo) - lo > 0) and (km.T(alpha, hi) - hi < 0) and (lo > q_lb) and (hi < q_ub)
        return bool(ok), lo, hi
    ok_l, ql_lo, ql_hi = tight(a_lb)
    ok_u, qu_lo, qu_hi = tight(a_ub)
    Gl = km.G(a_lb, box(ql_lo, ql_hi), km.R(a_lb, box(ql_lo, ql_hi)))
    Gu = km.G(a_ub, box(qu_lo, qu_hi), km.R(a_ub, box(qu_lo, qu_hi)))
    dGda = km.ElogPsi(Q)
    c3 = bool(ok_l and ok_u and (Gl > 0) and (Gu < 0) and (dGda < 0))
    Q0 = box(ql_lo, qu_hi)
    Psi0 = km.R(A, Q0)
    out = dict(kappa=kappa, alpha_lb=a_lb.str(30, radius=False), alpha_ub=a_ub.str(30, radius=False),
               q_lb=q_lb.str(30, radius=False), q_ub=q_ub.str(30, radius=False),
               sup_dT=str(dTmax), T_minus_q_at_lb=str(s_lb), T_minus_q_at_ub=str(s_ub),
               G_at_alpha_lb=str(Gl), G_at_alpha_ub=str(Gu), dG_dalpha=str(dGda),
               q0=str(Q0), psi0=str(Psi0), q0_lo=ql_lo.str(40, radius=False), q0_hi=qu_hi.str(40, radius=False),
               B1=bool(c1), B2=c2, B3=c3, seconds=round(time.time()-t0, 1))
    if verbose: print(json.dumps(out, indent=1), flush=True)
    return dict(A=A, Q0=Q0, Psi0=Psi0, a_lb=a_lb, a_ub=a_ub, ok=bool(c1 and c2 and c3), report=out)

FLOATS = {0.05: (0.7810743282532171, 0.5720535360129915), 0.13: (0.7059325315101108, 0.5847018482132774),
          0.0: (0.8330785997395772, 0.5639490800845581)}
if __name__ == "__main__":
    for k in [float(x) for x in sys.argv[1:]] or [0.05]:
        fp_box(k, *FLOATS[k])
