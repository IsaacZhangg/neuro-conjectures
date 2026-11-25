"""Certify Huang's Conditions 3.1, 3.3, 3.4 at a given kappa with Arb ball arithmetic."""
import sys, json, time
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arbgauss import phi, Psi, Ecal, gauss_E, ub
from flint import arb, acb, ctx

L = 14

def box(lo, hi):
    lo, hi = arb(lo), arb(hi)
    return lo.union(hi)

class KM:
    def __init__(self, kappa):
        self.k = arb(kappa)
    # ---- building blocks
    def P(self, psi):
        r = psi.sqrt()
        return gauss_E(lambda z, an: (r*z).tanh()**2, L, (1.0,), even=True)
    def dP(self, psi):  # P'(psi) = E[Z th sech^2]/sqrt(psi), integrand >= 0
        r = psi.sqrt()
        I = gauss_E(lambda z, an: z*(r*z).tanh()*(r*z).sech()**2, L, (0.0, 1.0), even=True)
        return I/r
    def _ua(self, q):
        s = (1-q).sqrt()
        a0 = ub((abs(self.k))/s); a1 = ub(q.sqrt()/s)
        return s, a0, a1
    def R(self, alpha, q):
        s, a0, a1 = self._ua(q)
        k, rq = self.k, q.sqrt()
        f = lambda z, an: Ecal((k - rq*z)/s)**2
        # E(u)^2 <= (|u|+1)^2 <= 2(a0+1)^2 + 2 a1^2 z^2
        I = gauss_E(f, L, (2*(a0+1)**2, 0.0, 2*a1**2))
        return alpha*I/(1-q)
    def dR(self, alpha, q):
        s, a0, a1 = self._ua(q)
        k, rq = self.k, q.sqrt()
        def f(z, an):
            u = (k - rq*z)/s
            E = Ecal(u); Ep = E*(E-u)
            return E*E/(1-q)**2 + E*Ep*(-z/(rq*s) + u/(1-q))/(1-q)
        # |E| <= |u|+1, 0<E'<1 on reals; crude polynomial bound
        C = ub(((a0+1)**2/(1-q)**2 + (a0+1+a1)*(1/(rq*s) + (a0+a1)/(1-q))/(1-q))*4)
        I = gauss_E(f, L, (C, C, C))
        return alpha*I
    def T(self, alpha, q):
        return self.P(self.R(alpha, q))
    def dT(self, alpha, q):
        return self.dP(self.R(alpha, q))*self.dR(alpha, q)
    def Elog2ch(self, psi):
        r = psi.sqrt()
        f = lambda z, an: r*z + (1 + (-2*r*z).exp()).log(analytic=an)
        return gauss_E(f, L, (1.0, ub(r)), even=True)
    def ElogPsi(self, q):
        s, a0, a1 = self._ua(q)
        k, rq = self.k, q.sqrt()
        f = lambda z, an: Psi((k - rq*z)/s).log(analytic=an)
        return gauss_E(f, L, (a0*a0/2 + a0 + 2, a0*a1 + a1, a1*a1/2))
    def G(self, alpha, q, psi):
        return -(1-q)*psi/2 + self.Elog2ch(psi) + alpha*self.ElogPsi(q)
    # ---- Condition 3.3 / 3.4 pieces
    def Esech4(self, psi):
        r = psi.sqrt()
        return gauss_E(lambda z, an: (r*z).sech()**4, L, (1.0,), even=True)
    def EFp(self, q, power):
        s, a0, a1 = self._ua(q)
        k, rq = self.k, q.sqrt()
        def f(z, an):
            u = (k - rq*z)/s
            E = Ecal(u)
            return (-(E*(E-u))/(1-q))**power
        return gauss_E(f, L, (ub(1/(1-q)**power),))
    def m(self, z0, psi):
        r = psi.sqrt()
        return gauss_E(lambda z, an: 1/(z0 + (r*z).cosh()**2), L, (ub(1/(1+z0)),), even=True)
    def Eg(self, q, mz):
        s, a0, a1 = self._ua(q)
        k, rq = self.k, q.sqrt()
        def f(z, an):
            u = (k - rq*z)/s
            E = Ecal(u); Ep = E*(E-u)
            return Ep/((1-q)*(1-Ep) + mz*Ep)
        mlo = float(mz.mid()) - float(mz.rad())
        assert mlo > 0
        return gauss_E(f, L, (1.0/mlo*1.0001,))
    def lam(self, z0, alpha, q, psi):
        mz = self.m(z0, psi)
        d0 = alpha*self.EFp(q, 1)
        return z0 - alpha*self.Eg(q, mz) - d0

def fixed_point_mid(km, alpha, q):
    for _ in range(400):
        qn = km.T(alpha, arb(q.mid()))
        qn = arb(qn.mid())
        if abs(float(qn - q)) < 1e-30: break
        q = qn
    return q

def certify(kappa, alpha_f, q_f, dalpha=1e-9, dq=1e-7, prec=128, zgrid=None, verbose=True):
    ctx.prec = prec
    km = KM(kappa)
    out = {"kappa": kappa}
    t0 = time.time()
    a_lb = arb(alpha_f) - arb(dalpha); a_ub = arb(alpha_f) + arb(dalpha)
    a_lb = arb(a_lb.mid()); a_ub = arb(a_ub.mid())
    A = box(a_lb, a_ub)
    q_lb = arb(arb(q_f - dq).mid()); q_ub = arb(arb(q_f + dq).mid())
    Q = box(q_lb, q_ub)
    # (i) contraction on the whole box
    dTmax = km.dT(A, Q)
    c1 = ub(dTmax) < 1
    # (ii) endpoint signs (T increasing in alpha)
    s_lb = km.T(a_lb, q_lb) - q_lb
    s_ub = km.T(a_ub, q_ub) - q_ub
    c2 = (s_lb > 0) and (s_ub < 0)
    # (iii) G* at alpha endpoints, using tight fixed-point enclosures
    def tight(alpha):
        qm = fixed_point_mid(km, alpha, arb(q_f))
        e = arb('1e-28')
        qm_lo, qm_hi = qm - e, qm + e
        qm_lo = arb(qm_lo.mid()); qm_hi = arb(qm_hi.mid())
        ok = (km.T(alpha, qm_lo) - qm_lo > 0) and (km.T(alpha, qm_hi) - qm_hi < 0) and (qm_lo > q_lb) and (qm_hi < q_ub)
        return ok, qm_lo, qm_hi
    ok_l, ql_lo, ql_hi = tight(a_lb)
    ok_u, qu_lo, qu_hi = tight(a_ub)
    Gl = km.G(a_lb, box(ql_lo, ql_hi), km.R(a_lb, box(ql_lo, ql_hi)))
    Gu = km.G(a_ub, box(qu_lo, qu_hi), km.R(a_ub, box(qu_lo, qu_hi)))
    dGda = km.ElogPsi(Q)   # dG*/dalpha = E log Psi(u) at q*, enclosed over Q
    c3 = ok_l and ok_u and (Gl > 0) and (Gu < 0) and (dGda < 0)
    Q0 = box(ql_lo, qu_hi)              # q* increasing in alpha
    Psi0 = km.R(A, Q0)
    out.update(alpha_lb=str(a_lb.mid()), alpha_ub=str(a_ub.mid()), q_lb=str(q_lb.mid()), q_ub=str(q_ub.mid()),
               sup_dT=str(dTmax), T_minus_q_at_lb=str(s_lb), T_minus_q_at_ub=str(s_ub),
               G_at_alpha_lb=str(Gl), G_at_alpha_ub=str(Gu), dG_dalpha=str(dGda),
               q0=str(Q0), psi0=str(Psi0), cond31=bool(c1 and c2 and c3))
    # Condition 3.3
    c33v = A*km.Esech4(Psi0)*km.EFp(Q0, 2)
    out.update(cond33_value=str(c33v), cond33=bool(c33v < 1))
    # Condition 3.4: choose zhat by a coarse midpoint scan, then certify at zhat
    am, qm, pm = arb(A.mid()), arb(Q0.mid()), arb(Psi0.mid())
    ctx.prec = 64
    best = None
    grid = zgrid or [-0.99 + 0.01*i for i in range(0, 300)]
    for z in grid:
        v = float(km.lam(arb(z), am, qm, pm).mid())
        if best is None or v < best[1]: best = (z, v)
    zf = best[0]
    for h in [0.004, 0.001, 0.0002]:
        cands = [zf + h*j for j in range(-5, 6) if zf + h*j > -1]
        vals = [(z, float(km.lam(arb(z), am, qm, pm).mid())) for z in cands]
        zf = min(vals, key=lambda t: t[1])[0]
    ctx.prec = prec
    zhat = arb(round(zf, 6))
    lamv = km.lam(zhat, A, Q0, Psi0)
    out.update(zhat=str(zhat), lambda_at_zhat=str(lamv), cond34=bool(lamv < 0))
    out["seconds"] = round(time.time()-t0, 1)
    if verbose: print(json.dumps(out, indent=1))
    return out

if __name__ == "__main__":
    vals = {0:(0.8330785997,0.5639490804),0.1:(0.7330166680,0.5800041202),0.25:(0.6092595657,0.6029620657),
            0.5:(0.4556645646,0.6384133980),0.75:(0.3478935315,0.6706271732),1.0:(0.2707000933,0.6998884161)}
    ks = [float(x) for x in sys.argv[1:]] or [0]
    for k in ks:
        certify(k, *vals[k])
