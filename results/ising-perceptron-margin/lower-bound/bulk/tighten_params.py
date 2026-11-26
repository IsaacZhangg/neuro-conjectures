# AGENT-WRITTEN (Claude), NOT INDEPENDENTLY REVIEWED.
# Tighter enclosure of psi* over the certified alpha-box.
# certify_1d gives psi0 = R(alpha-box, Q0) (width ~4e-7, from evaluating R on product boxes).
# Here: q*(alpha) and psi*(alpha) = P^{-1}(q*(alpha)) are both increasing in alpha on the box
# (T = P(R) is increasing in alpha, sup dT/dq < 1 on the certified box, and P is increasing),
# so psi*(alpha) lies in [psi*(alpha_lb), psi*(alpha_ub)], and each endpoint value is enclosed
# by R(alpha_end, tight q*(alpha_end)), with the tight q*(alpha_end) re-certified by sign changes
# of T(alpha_end, q) - q exactly as in certify_1d.tight().
import sys, json
sys.path.insert(0, '/home/claude/work')
from certify_1d import KM, fixed_point_mid, box
from arbgauss import ub
from flint import arb, ctx
ctx.prec = 128
pf = sys.argv[1]
P = json.load(open(pf))
km = KM(P['kappa'])
a_lb, a_ub = arb(P['alpha_lb']).mid(), arb(P['alpha_ub']).mid()
a_lb, a_ub = arb(a_lb), arb(a_ub)
q_lb, q_ub = arb(P['q_lb']), arb(P['q_ub'])
qf = arb(arb(P['q0']).mid())
def tight(alpha):
    qm = fixed_point_mid(km, alpha, qf)
    e = arb('1e-28')
    lo, hi = arb((qm - e).mid()), arb((qm + e).mid())
    ok = (km.T(alpha, lo) - lo > 0) and (km.T(alpha, hi) - hi < 0) and (lo > q_lb) and (hi < q_ub)
    return ok, lo, hi
okl, ql0, ql1 = tight(a_lb)
oku, qu0, qu1 = tight(a_ub)
psil = km.R(a_lb, box(ql0, ql1))
psiu = km.R(a_ub, box(qu0, qu1))
psi0 = box(psil, psiu)
q0 = box(ql0, qu1)
assert okl and oku
assert arb(P['psi0']).contains(psi0) or arb(P['psi0']).overlaps(psi0)
P2 = dict(P)
P2.update(q0=str(q0), psi0=str(psi0), psi_at_alpha_lb=str(psil), psi_at_alpha_ub=str(psiu),
          tight_q_alpha_lb=[str(ql0), str(ql1)], tight_q_alpha_ub=[str(qu0), str(qu1)],
          psi0_certify_1d=P['psi0'], q0_certify_1d=P['q0'],
          note='psi0 tightened by monotonicity of psi*(alpha); see tighten_params.py')
json.dump(P2, open(pf.replace('.json', '_tight.json'), 'w'), indent=1)
print(P['q0'], '->', q0); print(P['psi0'], '->', psi0)
