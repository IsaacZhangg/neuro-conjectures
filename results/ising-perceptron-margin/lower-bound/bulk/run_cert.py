# Agent-written (Claude), not independently reviewed.
# Certify the RS fixed point (alpha*, q*, psi*) at margin kappa with certify_1d (Arb ball arithmetic).
# kappa is passed as a decimal string so that arb(kappa) is an enclosure of the exact decimal.
# usage: run_cert.py kappa alpha_float q_float [dalpha dq tag]
import sys, json
sys.path.insert(0, '/home/claude/work')
from certify_1d import certify
k = sys.argv[1]
a = float(sys.argv[2]); q = float(sys.argv[3])
kw = {}
tag = ''
if len(sys.argv) > 5:
    kw = dict(dalpha=float(sys.argv[4]), dq=float(sys.argv[5]))
    tag = sys.argv[6] if len(sys.argv) > 6 else '_narrow'
out = certify(k, a, q, **kw)
out.update(kw)
json.dump(out, open(f'params_{k}{tag}.json', 'w'), indent=1)
