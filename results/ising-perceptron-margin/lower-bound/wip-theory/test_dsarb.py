# consistency test of dsarb.py against the float evaluator (NOT a certificate)
import sys, math, time
sys.path.insert(0, '/tmp/claude-0/-home-user/fa4b6461-5216-5a3d-ae79-fb39767abede/scratchpad/dingsun')
sys.path.insert(0, '/tmp/claude-0/-home-user/fa4b6461-5216-5a3d-ae79-fb39767abede/scratchpad/ds_theory')
import numpy as np
import ds_S as d
from dsarb import *
ctx.prec = 64
k = 0.05
ds = d.DS(k)
bx = Box(k, ds.alpha, ds.q, ds.psi)
t=time.time(); HM = h_consts(); print('H0,M1,M2', HM, time.time()-t)
for lam in [-0.3, 0.2, 0.9]:
    A = ds.A_of_lam(lam); a = math.log(A)
    t=time.time(); l = bx.ell(a); t1=time.time()-t
    t=time.time(); Lv = bx.Lfun(a); t2=time.time()-t
    t=time.time(); en = bx.entropy(a); t3=time.time()-t
    Hs = bx.Hstar()
    print('lam',lam,'ell',l,'L',Lv,'H via ent', -2*Hs+en, 'float H', ds.Hfun(lam), 'times',round(t1,3),round(t2,3),round(t3,3))
    xi = ds.xi; E = ds.Exi; w = d.ZW; c0 = math.sqrt(ds.psi*(1-ds.q)); a_ = ds.alpha
    for s in [0.0, 0.3]:
        c = math.sqrt((1-lam)/(1+lam)); th = E*s/c0
        ft = a_*float(np.sum(w*E*d.Emills(c*xi-th)))/math.sqrt(1-lam*lam)
        gb = (xi-lam*E)/math.sqrt(1-lam*lam)-th; Eg = d.Emills(gb)
        wj = float(np.sum(w*(Eg*Eg+th*Eg)))
        DQf = -ds.psi*(1-ds.q)/(1+lam)**2 + c0*s/((1+lam)*math.sqrt(1-lam*lam)) + ft - a_*lam/(1-lam*lam)*wj
        t=time.time(); DQa = bx.DP(lam, s); t4=time.time()-t
        print('   s',s,'DP arb',DQa,'float',DQf,'time',round(t4,3))
    # cell (ball) evaluation timing
    t=time.time(); DQb = bx.DP(hull(lam, lam+0.002), 0.3); print('   ball DP', DQb, round(time.time()-t,3))
print('Q0(0.3)', bx.Q0(0.3), 'float', 0.3**2/2 - c0*0.3 + a_*float(np.sum(w*(d.logPsi(xi-E*0.3/c0)-d.logPsi(xi)))))
print('lam_min', bx.lam_min(), ds.lam_min)
for lam in [0.9, 0.95]:
    t=time.time(); print('K', lam, bx.Knear1(lam, HM), bx.Knear1(lam, HM, improved=False), round(time.time()-t,3))
print('Kunif 0.05', bx.Kunif(0.05, HM))
print('F1prime(0)', bx.F1prime(0), 'L(0)', bx.Lfun(0))
