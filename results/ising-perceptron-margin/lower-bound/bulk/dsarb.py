# AGENT-WRITTEN (Claude), NOT INDEPENDENTLY REVIEWED.
#
# Rigorous (Arb ball arithmetic, python-flint) enclosures of the pieces of the
# Ding--Sun second-moment exponent S(lambda) (arXiv:1809.07742 v1, eqs (34)-(48)),
# at general margin kappa, for parameters (alpha, q, psi) given as Arb balls.
#
# Notation:  Psi = upper Gaussian tail, E = phi/Psi,  xi_z = (kappa - sqrt(q) z)/sqrt(1-q),
#            kk = sqrt(psi) sqrt(1-q),  f(lam) = (1-lam)/(1+lam),  w = sqrt(1-lam^2).
#
#   H(A)    = -2 Hstar + E Gamma(sqrt(psi) Z, D(A))          (entropy part, 1-D)
#   ell(A)  = E D_{sqrt(psi)Z}(A) / (1-q)                     (lambda as a function of A, 1-D)
#   H'(lam) = -(1-q) ln A(lam) / 2      (exact; H is concave in lam, see README)
#   PB(lam) = P(lam) + B(lam, s)
#           = -Pstar + psi(1-q) f/2 + s^2/2 - kk s sqrt(f) + I_s(lam)   (the I_0 terms cancel)
#   I_s(lam)= alpha E_z E_{nu ~ phi_xi} ln Psi(arg),  arg = (xi - lam nu)/w - E(xi) s / kk
#
# The 2-D expectation is a nested acb.integral: outer z in [-L, L] split into unit pieces,
# inner nu in [xi, c_piece] via nu = xi + (c - xi) t, t in [0, 1]; both truncations are
# bounded by explicit polynomial envelopes (functions tailbound_*).
from flint import arb, acb, ctx
import math
import time

ctx.prec = 128
SQ2 = arb(2).sqrt()
SQ2PI = (2 * arb.pi()).sqrt()
LN4 = arb(4).log()


def Psi(x):
    return (x / SQ2).erfc() / 2


INV_SQ2PI = 1 / SQ2PI


def _Ps_direct(x):
    return Psi(x) * (x * x / 2).exp()


def Ps(x):
    """scaled tail g(x) = Psi(x) exp(x^2/2).
    Arb's erfc on complex balls with large positive real part is very pessimistic, so when
    Re x >= c > 1/2 on the ball we evaluate g at the (exact) midpoint and add rad * sup|g'|:
    g(x) = (2pi)^{-1/2} int_0^inf exp(-x t - t^2/2) dt gives |g'(x)| <= (2pi)^{-1/2}/c^2 there.
    Elsewhere the direct ball evaluation is used."""
    if not isinstance(x, acb):
        x = acb(x)
    if x.is_exact():
        return _Ps_direct(x)
    c = x.real.lower()
    if not (c > 0.5):
        return _Ps_direct(x)
    m = x.mid()
    gm = _Ps_direct(m)
    r = arb(ub(x.rad()))
    D = INV_SQ2PI / (c * c)
    e = arb(0, ub(r * D))
    return gm + acb(e, e)


def lnPsi(x, an):
    if float(x.real.mid()) > 1.0:
        return -(x * x) / 2 + Ps(x).log(analytic=an)
    return Psi(x).log(analytic=an)


def Ecal(x):
    if float(acb(x).real.mid()) > 1.0:
        return 1 / (SQ2PI * Ps(x))
    return phi(x) / Psi(x)


def E_Ep(x):
    """(E(x), E'(x)) with E' = E (E - x).  For a non-exact ball with Re x >= c > 1/2, uses
    E = 1/(sqrt(2pi) g), E' = -g'/(sqrt(2pi) g^2), where g = Ps and g' = x g - 1/sqrt(2pi) are enclosed
    by midpoint values plus rad * sup|g'|, rad * sup|g''| with |g^(k)(x)| <= (2pi)^{-1/2} k!/c^(k+1)
    (from g(x) = (2pi)^{-1/2} int_0^inf e^{-xt - t^2/2} dt).  This avoids the cancellation in E - x."""
    if not isinstance(x, acb):
        x = acb(x)
    c = x.real.lower()
    if (not x.is_exact()) and c > 0.5:
        m = x.mid()
        gm = _Ps_direct(m)
        g1m = m * gm - INV_SQ2PI
        r = arb(ub(x.rad()))
        e1 = arb(0, ub(r * INV_SQ2PI / (c * c)))
        e2 = arb(0, ub(r * 2 * INV_SQ2PI / (c * c * c)))
        g = gm + acb(e1, e1)
        g1 = g1m + acb(e2, e2)
        E = 1 / (SQ2PI * g)
        return E, -g1 / (SQ2PI * g * g)
    E = Ecal(x)
    return E, E * (E - x)


def phi(x):
    return (-(x * x) / 2).exp() / SQ2PI


def ub(x):
    """float upper bound of a real arb (outward rounded)."""
    x = arb(x)
    if not x.is_finite():
        return float('inf')
    v = float(x.mid()) + float(x.rad())
    return v + abs(v) * 1e-14 + 1e-300


def lb(x):
    x = arb(x)
    if not x.is_finite():
        return float('-inf')
    v = float(x.mid()) - float(x.rad())
    return v - abs(v) * 1e-14 - 1e-300


def box(lo, hi):
    return arb(lo).union(arb(hi))


def T2(L):
    """int_{|z|>L} z^2 phi(z) dz = 2 (L phi(L) + Psi(L))  (exact identity)."""
    L = arb(L)
    return 2 * (L * phi(L) + Psi(L))


def T0(L):
    """int_{|z|>L} phi = 2 Psi(L)."""
    return 2 * Psi(arb(L))


class Model:
    def __init__(self, kappa, alpha, q, psi, L=10, Lh=10, tol=1e-14):
        self.k = arb(kappa)
        self.a = arb(alpha)
        self.q = arb(q)
        self.psi = arb(psi)
        self.L = L
        self.Lh = Lh
        self.tol = tol
        self.rq = self.q.sqrt()
        self.w1q = (1 - self.q).sqrt()
        self.rpsi = self.psi.sqrt()
        self.kk = self.rpsi * self.w1q
        self.xi0 = self.k / self.w1q          # xi = xi0 - xi1 z
        self.xi1 = self.rq / self.w1q
        assert self.k >= 0
        self._Hstar = None
        self.ball_eval_limit = 3000
        self.range_eval_limit = 600
        self.counts = [0, 0]
        self._Pstar = None

    # ------------------------------------------------------------------ 1-D helpers
    def _int_even(self, f, L, pieces, tail):
        """2 int_0^L f(z) phi(z) dz + tail (tail = rigorous bound for |z| > L part)."""
        g = lambda z, an: f(z, an) * phi(z)
        h = arb(L) / pieces
        I = acb(0)
        for i in range(pieces):
            I += acb.integral(g, h * i, h * (i + 1), rel_tol=self.tol, abs_tol=self.tol * 1e-3)
        return 2 * I.real + arb(0, ub(tail))

    def _int_full(self, f, L, pieces, tail):
        g = lambda z, an: f(z, an) * phi(z)
        h = arb(2 * L) / pieces
        I = acb(0)
        for i in range(pieces):
            I += acb.integral(g, -L + h * i, -L + h * (i + 1), rel_tol=self.tol, abs_tol=self.tol * 1e-3)
        return I.real + arb(0, ub(tail))

    def Hstar(self):
        if self._Hstar is None:
            r = self.rpsi
            # ln 2ch(x) = x + ln(1+e^{-2x}) for x >= 0;  |ln 2ch x| <= |x| + ln 2
            f = lambda z, an: r * z + (1 + (-2 * r * z).exp()).log(analytic=an)
            L = self.Lh
            # tail: int_{|z|>L} (r|z| + 1) phi <= (r L + 1)/L^2 * T2(L)  (r|z|+1 <= (rL+1) z^2/L^2 for |z|>=L>=1)
            tail = (r * L + 1) / (L * L) * T2(L)
            E = self._int_even(f, L, 10, tail)
            self._Hstar = -self.psi * (1 - self.q) + E
        return self._Hstar

    def Pstar(self):
        if self._Pstar is None:
            k, rq, s = self.k, self.rq, self.w1q
            f = lambda z, an: lnPsi((k - rq * z) / s, an)
            L = self.L
            X = self.xi0 + self.xi1 * L
            # |ln Psi(x)| <= x^2/2 + |x| + 2, |xi| <= xi0 + xi1 |z|
            tail = (X * X / 2 + X + 2) / (L * L) * T2(L)
            E = self._int_full(f, L, 20, tail)
            self._Pstar = self.psi * (1 - self.q) / 2 + self.a * E
        return self._Pstar

    def _hpieces(self, z, A, an, need_G=True):
        """m = tanh(sqrt(psi) z) for z >= 0 (even integrands);  returns D, Gamma.
        D = (A^2-1)(1-m^2)^2/(Delta+1)^2,  Delta = sqrt(m^2 + A^2(1-m^2)) = sqrt(1 + (A^2-1)(1-m^2)).
        Gamma is computed from the weights (A B, 1, 1, A/B) proportional to (p11, p12, p21, p22),
        A B = (Delta+m)/(1-m), A/B = A^2 (1-m)/(Delta+m)  (the Q_{A,B} form, eq. (127)-(129));
        this avoids the cancellations in p22 (A small) and p12 (A large)."""
        x = self.rpsi * z
        omm = 2 / ((2 * x).exp() + 1)        # 1 - m
        opm = 2 - omm                        # 1 + m
        m = 1 - omm
        c = omm * opm                        # 1 - m^2
        A2 = A * A
        Dl = (1 + (A2 - 1) * c).sqrt(analytic=an)
        D = (A2 - 1) * c * c / (Dl + 1) ** 2
        if not need_G:
            return D, None
        u1 = (Dl + m) / omm
        u4 = A2 * omm / (Dl + m)
        Zt = u1 + 2 + u4
        G = Zt.log(analytic=an) - (u1 * u1.log(analytic=an) + u4 * u4.log(analytic=an)) / Zt
        return D, G

    def ell(self, A):
        A = arb(A)
        f = lambda z, an: self._hpieces(z, A, an, need_G=False)[0]
        L = self.Lh
        tail = T0(L)                         # |D| <= 1 - m^2 <= 1
        return self._int_even(f, L, 10, tail) / (1 - self.q)

    def Hfun(self, A):
        A = arb(A)
        f = lambda z, an: self._hpieces(z, A, an)[1]
        L = self.Lh
        tail = LN4 * T0(L)                   # 0 <= Gamma <= ln 4
        return -2 * self.Hstar() + self._int_even(f, L, 10, tail)

    def lam_min(self):
        """ell(0) = -E (1-|m|)^2/(1-q)."""
        r = self.rpsi
        f = lambda z, an: -(2 / ((2 * r * z).exp() + 1)) ** 2
        return self._int_even(f, self.Lh, 10, T0(self.Lh)) / (1 - self.q)

    def Hprime(self, A):
        """dH/dlam at lam = ell(A)."""
        return -(1 - self.q) * arb(A).log() / 2

    # ------------------------------------------------------------------ 2-D part
    def _envelope(self, kind, X, lam_abs, w, s_abs):
        """coefficients e0, e1, e2 (arb, >= 0) with |g| <= e0 + e1 |nu| + e2 nu^2 whenever |xi| <= X."""
        EX = X + arb('0.8')
        aB = X / w + EX * s_abs / self.kk
        bB = lam_abs / w
        if kind == 0:
            return aB * aB + aB + 2, bB, bB * bB
        w3 = w ** 3
        if kind == 1:
            return (aB + arb('0.8')) * lam_abs * X / w3, ((aB + arb('0.8')) + bB * lam_abs * X) / w3, bB / w3
        w5 = w ** 5
        w6 = w ** 6
        p0 = X / w3 + 3 * lam_abs * lam_abs * X / w5
        p1 = 3 * lam_abs / w5
        e0 = 2 * lam_abs * lam_abs * X * X / w6 + (aB + arb('0.8')) * p0
        e1 = (aB + arb('0.8')) * p1 + bB * p0
        e2 = 2 / w6 + bB * p1
        return e0, e1, e2

    def Igen(self, kind, lam, s, pieces_per_unit=1, inner_tol=None, outer_tol=None, outer_abs=None, verbose=False, L=None):
        """alpha E_z E_{nu~phi_xi} g_kind, with g_0 = ln Psi(arg), g_1 = d/dlam g_0, g_2 = d^2/dlam^2 g_0.
        lam may be a (real) ball; the result then encloses the value for every lam in the ball."""
        lam = arb(lam)
        s = arb(s)
        w = (1 - lam * lam).sqrt()
        assert w > 0
        k, rq, s1q, kk = self.k, self.rq, self.w1q, self.kk
        w3 = w ** 3
        w5 = w ** 5
        itol = inner_tol or self.tol
        lam_c, w_c, s_c, kk_c = acb(lam), acb(w), acb(s), acb(kk)
        w3c, w5c = acb(w3), acb(w5)

        def make_outer(cc):
            ccc = acb(cc)

            def outer(z, an_out):
                self.counts[0] += 1
                xi = (k - rq * z) / s1q
                Pxi = Ps(xi)                        # Psi(xi) e^{xi^2/2}
                Exi = 1 / (SQ2PI * Pxi)
                span = ccc - xi
                a0 = xi / w_c - Exi * s_c / kk_c

                def inner(t, an_in):
                    self.counts[1] += 1
                    an = an_in or an_out
                    u = span * t
                    nu = xi + u
                    # phi(nu)/Psi(xi) = exp(-(nu^2 - xi^2)/2) / (sqrt(2 pi) Ps(xi))
                    wt = (-(u * (2 * xi + u)) / 2).exp()
                    x = a0 - lam_c * nu / w_c
                    if kind == 0:
                        g = lnPsi(x, an)
                    elif kind == 'p01':
                        # pack g_0 + i g_1 (both real on the real axis): one nested integral gives both
                        g = lnPsi(x, an) + acb(0, 1) * (-Ecal(x) * (lam_c * xi - nu) / w3c)
                    else:
                        xl = (lam_c * xi - nu) / w3c
                        if kind == 1:
                            g = -Ecal(x) * xl
                        else:
                            Ex, Epx = E_Ep(x)
                            xll = xi / w3c + 3 * lam_c * (lam_c * xi - nu) / w5c
                            g = -Epx * xl * xl - Ex * xll
                    return wt * g

                if ub(z.rad()) < 1e-12:
                    J = acb.integral(inner, 0, 1, rel_tol=itol, abs_tol=itol * 1e-3)
                else:
                    # z is a wide (complex) ball: the outer integrator only needs an enclosure here, for
                    # its error bound.  Use a range enclosure: int_0^1 h(t) dt is contained in
                    # (1/n) sum_j h([j/n, (j+1)/n]) (h evaluated on real sub-balls of [0, 1]).
                    J = None
                    for nsub in (8, 64):
                        acc = acb(0)
                        for j in range(nsub):
                            tb = acb(arb((2 * j + 1), 1) / (2 * nsub))
                            acc += inner(tb, False)
                        acc = acc / nsub
                        if acc.is_finite():
                            J = acc
                            break
                    if J is None:
                        J = acb(arb('nan'))
                return J * span / (SQ2PI * Pxi) * phi(z)
            return outer

        L = L or self.L        # truncation |z| <= L (integer); |z| > L is bounded by the envelope
        n = 2 * L * pieces_per_unit
        h = arb(2 * L) / n
        tot = acb(0)
        envk = [0, 1] if kind == 'p01' else [kind]
        nutail = [arb(0) for _ in envk]
        lam_abs = arb(ub(abs(lam)))
        w_lo = arb(lb(w))
        s_abs = arb(ub(abs(s)))
        for i in range(n):
            za = -L + h * i
            zb = za + h
            xi_a = (k - rq * za) / s1q          # sup of xi on the piece
            xi_b = (k - rq * zb) / s1q
            cc = math.ceil(max(ub(xi_a), 0.0)) + 10
            _t = time.time(); _c = list(self.counts)
            tot += acb.integral(make_outer(cc), za, zb, rel_tol=(outer_tol or self.tol),
                                abs_tol=(outer_abs or self.tol * 1e-3))
            if verbose:
                print('piece', i, float(za.mid()), round(time.time() - _t, 2), self.counts[0] - _c[0], self.counts[1] - _c[1], tot, flush=True)
            # inner truncation nu > cc on this piece
            X = arb(max(ub(abs(xi_a)), ub(abs(xi_b))))
            C = arb(cc)
            mass = Psi(za) - Psi(zb)
            for j, kd in enumerate(envk):
                e0, e1, e2 = self._envelope(kd, X, lam_abs, w_lo, s_abs)
                num = e0 * Psi(C) + e1 * phi(C) + e2 * (C * phi(C) + Psi(C))
                nutail[j] += num / Psi(arb(ub(xi_a))) * mass   # Psi(xi) >= Psi(sup xi)
        # outer truncation |z| > L
        X = self.xi0 + self.xi1 * L
        res = []
        for j, kd in enumerate(envk):
            e0, e1, e2 = self._envelope(kd, X, lam_abs, w_lo, s_abs)
            PL = e0 + e1 * (X + arb('1.6')) + e2 * (1 + X * X + arb('0.8') * X)
            ztail = PL / (L * L) * T2(L)
            if verbose:
                print('kind', kd, 'nutail', nutail[j], 'ztail', ztail, 'core', tot)
            t = ub(nutail[j]) + ub(ztail)
            part = tot.real if j == 0 else tot.imag
            res.append(self.a * (part + arb(0, t)))
        return res[0] if len(res) == 1 else tuple(res)

    def I2_range(self, lam, s, hz=0.1, L=6, inner_tol=1e-3, verbose=False):
        """Enclosure of I_s''(lam) = alpha E_z E_nu g_2 valid for every lam in the (real) ball lam.
        Outer z-integral by a range (Riemann) enclosure: int_{za}^{zb} F(z) dz is contained in
        (zb - za) F([za, zb]), where F([za, zb]) is the inner acb.integral evaluated with z a real
        ball (acb.integral with ball parameters encloses the integral for every parameter value in
        the ball).  This avoids evaluating the inner integral at complex z (the outer integrator's
        error bound), which is where the nested ball computation blows up.  Truncations as in Igen."""
        lam = arb(lam)
        s = arb(s)
        w = (1 - lam * lam).sqrt()
        assert w > 0
        k, rq, s1q, kk = self.k, self.rq, self.w1q, self.kk
        lam_c, w_c, s_c, kk_c = acb(lam), acb(w), acb(s), acb(kk)
        w3c, w5c = acb(w ** 3), acb(w ** 5)
        n = int(round(2 * L / hz))
        h = arb(2 * L) / n
        tot = arb(0)
        nutail = arb(0)
        lam_abs = arb(ub(abs(lam)))
        w_lo = arb(lb(w))
        s_abs = arb(ub(abs(s)))
        t0 = time.time()
        for i in range(n):
            za = -L + h * i
            zb = za + h
            zr = box(za, zb)
            z = acb(zr)
            xi_a = (k - rq * za) / s1q
            xi_b = (k - rq * zb) / s1q
            cc = acb(math.ceil(max(ub(xi_a), 0.0)) + 10)
            xi = (k - rq * z) / s1q
            Pxi = Ps(xi)
            Exi = 1 / (SQ2PI * Pxi)
            span = cc - xi
            a0 = xi / w_c - Exi * s_c / kk_c

            def inner(t, an):
                self.counts[1] += 1
                u = span * t
                nu = xi + u
                wt = (-(u * (2 * xi + u)) / 2).exp()
                x = a0 - lam_c * nu / w_c
                xl = (lam_c * xi - nu) / w3c
                Ex, Epx = E_Ep(x)
                xll = xi / w3c + 3 * lam_c * (lam_c * xi - nu) / w5c
                return wt * (-Epx * xl * xl - Ex * xll)

            J = acb.integral(inner, 0, 1, rel_tol=inner_tol, abs_tol=inner_tol * 1e-2,
                             eval_limit=self.range_eval_limit)
            if not J.is_finite():
                # acb.integral gave up (ball widths dominate): range enclosure on real t-sub-balls,
                # int_0^1 h dt is contained in (1/N) sum_j h([j/N, (j+1)/N]).
                for N in (64, 512):
                    acc = acb(0)
                    for j in range(N):
                        acc += inner(acb(arb(2 * j + 1, 1) / (2 * N)), False)
                    if acc.is_finite():
                        J = acc / N
                        self.counts[0] += 1
                        break
            F = (J * span / (SQ2PI * Pxi) * phi(z)).real
            if not F.is_finite():
                return arb('nan')
            tot += h * F
            X = arb(max(ub(abs(xi_a)), ub(abs(xi_b))))
            C = arb(math.ceil(max(ub(xi_a), 0.0)) + 10)
            mass = Psi(za) - Psi(zb)
            e0, e1, e2 = self._envelope(2, X, lam_abs, w_lo, s_abs)
            num = e0 * Psi(C) + e1 * phi(C) + e2 * (C * phi(C) + Psi(C))
            nutail += num / Psi(arb(ub(xi_a))) * mass
        X = self.xi0 + self.xi1 * L
        e0, e1, e2 = self._envelope(2, X, lam_abs, w_lo, s_abs)
        PL = e0 + e1 * (X + arb('1.6')) + e2 * (1 + X * X + arb('0.8') * X)
        ztail = PL / (L * L) * T2(L)
        if verbose:
            print('I2_range', lam, 'core', tot, 'nutail', nutail, 'ztail', ztail, round(time.time() - t0, 1), flush=True)
        return self.a * (tot + arb(0, ub(nutail) + ub(ztail)))

    def I2_real(self, lam, s, hz=0.05, L=6, tau=0.01, hmin=1e-3, umax=9, verbose=False):
        """Enclosure of I_s''(lam) valid for every lam in the real ball lam, by a 2-D range
        (interval Riemann) enclosure in REAL ball arithmetic.  Uses nu = xi + u, u >= 0, with the
        conditional density phi(nu)/Psi(xi) = exp(-xi u - u^2/2)/(sqrt(2 pi) Ps(xi)), and on every cell
        int_cell f is contained in |cell| * f(cell).  The z-cells use the exact Gaussian mass
        Psi(za) - Psi(zb) (phi > 0, mean value theorem).  E = phi/Psi and E' on a real ball [a, b]
        come from monotonicity only (classical: E increasing, E - x decreasing and positive, i.e.
        0 < E' < 1, Sampford 1953):  E in [E(a), E(b)],  E - x in [E(b) - b, E(a) - a],
        E' = E (E - x) in the product.  Ps = Psi e^{x^2/2} is decreasing (Ps' = Ps (x - E) < 0)."""
        lam = arb(lam)
        s = arb(s)
        w = (1 - lam * lam).sqrt()
        assert w > 0
        k, rq, s1q, kk = self.k, self.rq, self.w1q, self.kk
        w3 = w ** 3
        w5 = w ** 5
        oml = 1 - lam
        n = int(round(2 * L / hz))
        h = arb(2 * L) / n
        lam_abs = arb(ub(abs(lam)))
        w_lo = arb(lb(w))
        s_abs = arb(ub(abs(s)))
        cacheE = {}

        def Ept(a):
            v = cacheE.get(a)
            if v is None:
                v = Ecal(acb(arb(a))).real
                cacheE[a] = v
            return v

        def E_rng(x):
            a, b = lb(x), ub(x)
            if not (math.isfinite(a) and math.isfinite(b)):
                return None
            Ea, Eb = Ept(a), Ept(b)
            E = Ea.union(Eb)
            D = (Eb - arb(b)).union(Ea - arb(a))
            return E, E * D

        tot = arb(0)
        nutail = arb(0)
        nevals = 0
        t0 = time.time()
        for i in range(n):
            za = -L + h * i
            zb = za + h
            xa = (k - rq * za) / s1q          # xi at za (the larger)
            xb = (k - rq * zb) / s1q
            xi = xa.union(xb)
            xlo, xhi = lb(xi), ub(xi)
            # Ps(xi) decreasing: Ps(xi) in [Ps(xhi), Ps(xlo)];  E(xi) = 1/(sqrt(2pi) Ps)
            Ps_r = (Ps(acb(arb(xhi))).real).union(Ps(acb(arb(xlo))).real)
            Exi = Ept(xlo).union(Ept(xhi))
            a0 = xi * oml / w - Exi * s / kk           # x = a0 - lam u / w
            U = math.ceil(max(0.0, -xlo)) + umax
            mass = Psi(za) - Psi(zb)

            def cell(u0, u1):
                u = arb(u0).union(arb(u1))
                wt = (-(u * (xi + u / 2))).exp()
                x = a0 - lam * u / w
                r = E_rng(x)
                if r is None:
                    return arb('nan')
                Ex, Epx = r
                dd = xi * oml + u                        # nu - lam xi
                xl2 = (dd * dd) / (w3 * w3)
                xll = xi / w3 - 3 * lam * dd / w5
                return (u1 - u0) * wt * (-Epx * xl2 - Ex * xll)

            # adaptive range sum over u in [0, U]
            acc = arb(0)
            stack = [(0.0, float(U), None)]
            while stack:
                u0, u1, v = stack.pop()
                if v is None:
                    v = cell(u0, u1)
                    nevals += 1
                hw = u1 - u0
                if (not v.is_finite()) or (ub(v.rad()) > tau * hw and hw > hmin):
                    um = 0.5 * (u0 + u1)
                    v1 = cell(u0, um); v2 = cell(um, u1)
                    nevals += 2
                    if v.is_finite() and v1.is_finite() and v2.is_finite() and \
                            ub(v1.rad()) + ub(v2.rad()) > 0.9 * ub(v.rad()) and hw < 0.05:
                        acc += v          # parameter-width dominated: refinement does not help
                        continue
                    stack.append((u0, um, v1)); stack.append((um, u1, v2))
                    continue
                acc += v
            F = acc / (SQ2PI * Ps_r)
            if not F.is_finite():
                return arb('nan')
            tot += mass * F
            # tail nu > xi + U (>= xlo + U >= 0 ... ) : envelope
            X = arb(max(abs(xlo), abs(xhi)))
            C = arb(xlo + U)
            e0, e1, e2 = self._envelope(2, X, lam_abs, w_lo, s_abs)
            num = e0 * Psi(C) + e1 * phi(C) + e2 * (C * phi(C) + Psi(C))
            nutail += num / Psi(arb(xhi)) * mass
        X = self.xi0 + self.xi1 * L
        e0, e1, e2 = self._envelope(2, X, lam_abs, w_lo, s_abs)
        PL = e0 + e1 * (X + arb('1.6')) + e2 * (1 + X * X + arb('0.8') * X)
        ztail = PL / (L * L) * T2(L)
        if verbose:
            print('I2_real', lam, 'core', tot, 'nutail', nutail, 'ztail', ztail, 'evals', nevals,
                  round(time.time() - t0, 1), flush=True)
        return self.a * (tot + arb(0, ub(nutail) + ub(ztail)))

    def I2_y(self, lam, s, hz=0.05, L=6, tau=0.01, hmin=1e-3, verbose=False):
        """Enclosure of I_s''(lam) for every lam in the real ball lam (needs lam > 0), using the
        boundary-layer variable nu = xi + w y (y >= 0), in which
            E_{nu ~ phi_xi} ln Psi(arg) = E(xi) int_0^inf K dy,
            K = w e^Q h(x),  Q = -xi w y - w^2 y^2/2,  x = xi r - lam y - c,
        with h = ln Psi, r = sqrt((1-lam)/(1+lam)), c = E(xi) s/kk, w = sqrt(1-lam^2).  The y-limits do
        not depend on lam, so I'' = alpha E_z E(xi) int_0^inf K'' dy with
            K'' = e^Q [ (w'' + 2 w'Q' + w Q'' + w Q'^2) h + 2 (w' + w Q') x' h' + w (x'^2 h'' + x'' h') ].
        This form has no catastrophic cancellation as lam -> 1.  h, h' = -E, h'' = -E' on a real ball
        from monotonicity (as in I2_real).  Range (interval Riemann) enclosure over z-cells (exact
        Gaussian mass) and adaptive y-cells on [0, Y]; the tail y > Y is bounded using
        |h|, |h'|, |h''| <= 2 (t+1) phi(t) for x = -t, t >= 1, and a ratio (geometric) bound."""
        lam = arb(lam)
        s = arb(s)
        assert lam > 0
        w = (1 - lam * lam).sqrt()
        assert w > 0
        k, rq, s1q, kk = self.k, self.rq, self.w1q, self.kk
        w1 = -lam / w
        w2 = -1 / w ** 3
        f = (1 - lam) / (1 + lam)
        r = f.sqrt()
        f1 = -2 / (1 + lam) ** 2
        f2 = 4 / (1 + lam) ** 3
        r1 = f1 / (2 * r)
        r2 = f2 / (2 * r) - f1 * f1 / (4 * r ** 3)
        lam_lo = arb(lb(lam))
        scale = 1 + ub(1 / w ** 3)
        n = int(round(2 * L / hz))
        h = arb(2 * L) / n
        lam_abs = arb(ub(abs(lam)))
        w_lo = arb(lb(w))
        s_abs = arb(ub(abs(s)))
        cache = {}

        def hpt(a):
            v = cache.get(a)
            if v is None:
                xa = acb(arb(a))
                v = (lnPsi(xa, False).real, Ecal(xa).real)
                cache[a] = v
            return v

        def h_rng(x):
            a, b = lb(x), ub(x)
            if not (math.isfinite(a) and math.isfinite(b)):
                return None
            (ha, Ea), (hb, Eb) = hpt(a), hpt(b)
            E = Ea.union(Eb)
            D = (Eb - arb(b)).union(Ea - arb(a))
            return ha.union(hb), -E, -(E * D)

        def Kpp(xi, c, y):
            Q = -(xi * w * y) - (w * y) ** 2 / 2
            Q1 = -(xi * w1 * y) - w * w1 * y * y
            Q2 = -(xi * w2 * y) - (w1 * w1 + w * w2) * y * y
            x = xi * r - lam * y - c
            x1 = xi * r1 - y
            x2 = xi * r2
            hr = h_rng(x)
            if hr is None:
                return arb('nan')
            h0, h1, h2 = hr
            return Q.exp() * ((w2 + 2 * w1 * Q1 + w * Q2 + w * Q1 * Q1) * h0
                              + 2 * (w1 + w * Q1) * x1 * h1 + w * (x1 * x1 * h2 + x2 * h1))

        tot = arb(0)
        ytail = arb(0)
        nevals = 0
        t0 = time.time()
        for i in range(n):
            za = -L + h * i
            zb = za + h
            xa = (k - rq * za) / s1q
            xb = (k - rq * zb) / s1q
            xi = xa.union(xb)
            xlo, xhi = lb(xi), ub(xi)
            Exi = hpt(xlo)[1].union(hpt(xhi)[1])
            c = Exi * s / kk
            mass = Psi(za) - Psi(zb)
            # truncation point Y:  for y >= Y,  t = -x >= lam y + b >= T,  b = lower bound of (c - xi r)
            a_ = max(0.0, -xlo) * ub(w)
            b_ = lb(c - xi * r)
            T = max(9.0, (2.0 * a_ + 2.0) / float(lam_lo) + 9.0)
            Y = max(1.0, math.ceil((T - b_) / float(lam_lo)))
            tolc = tau * scale / ub(Exi)
            acc = arb(0)
            stack = [(0.0, float(Y), None)]
            while stack:
                y0, y1, v = stack.pop()
                if v is None:
                    v = (y1 - y0) * Kpp(xi, c, arb(y0).union(arb(y1)))
                    nevals += 1
                hw = y1 - y0
                if (not v.is_finite()) or (ub(v.rad()) > tolc * hw and hw > hmin):
                    ym = 0.5 * (y0 + y1)
                    v1 = (ym - y0) * Kpp(xi, c, arb(y0).union(arb(ym)))
                    v2 = (y1 - ym) * Kpp(xi, c, arb(ym).union(arb(y1)))
                    nevals += 2
                    if v.is_finite() and v1.is_finite() and v2.is_finite() and \
                            ub(v1.rad()) + ub(v2.rad()) > 0.9 * ub(v.rad()) and hw < 0.05:
                        acc += v
                        continue
                    stack.append((y0, ym, v1)); stack.append((ym, y1, v2))
                    continue
                acc += v
            if not acc.is_finite():
                return arb('nan')
            tot += mass * Exi * acc
            # tail y > Y:  |K''| <= F(y) := e^{a y} P(y) 2 (t+1) phi(t),  t = lam_lo y + b  (>= T >= 3 > 1),
            # with P(y) = |A0| + |A1| + |A2| evaluated by ball arithmetic.  For y >= Y,
            # F(y+1)/F(y) <= rho := (1 + 1/(1+Y))^5 * (1 + lam/(T+1)) * e^{a - lam^2/2 - lam T}
            # (P has degree <= 4 in y with |coefficients| absorbed by (1+y)^4 up to a constant: we bound
            #  P(y) <= P4 (1+y)^4 with P4 = sup over y >= 0 of P(y)/(1+y)^4, see below).
            X = arb(max(abs(xlo), abs(xhi)))
            ax = X
            # coefficient bounds (all >= 0):  |Q1| <= ax|w1| y + w|w1| y^2,  |Q2| <= ax|w2| y + (w1^2 + w|w2|) y^2,
            # |x1| <= ax|r1| + y,  |x2| <= ax|r2|.  Each is <= (const) (1+y)^k; collect P(y) <= P4 (1+y)^4.
            aw1, aw2, ar1, ar2 = abs(w1), abs(w2), abs(r1), abs(r2)
            cQ1 = ax * aw1 + w * aw1                  # |Q1| <= cQ1 (1+y)^2
            cQ2 = ax * aw2 + aw1 * aw1 + w * aw2      # |Q2| <= cQ2 (1+y)^2
            cx1 = ax * ar1 + 1                        # |x1| <= cx1 (1+y)
            cx2 = ax * ar2
            P4 = (aw2 + 2 * aw1 * cQ1 + w * cQ2 + w * cQ1 * cQ1) + 2 * (aw1 + w * cQ1) * cx1 + w * (cx1 * cx1 + cx2)
            Yb = arb(Y)
            Tb = lam_lo * Yb + arb(b_)
            if not (Tb >= T - 1e-9):
                return arb('nan')
            aa = arb(a_)
            rho = (1 + 1 / (1 + Yb)) ** 5 * (1 + lam_abs / (Tb + 1)) * (aa - lam_lo * lam_lo / 2 - lam_lo * Tb).exp()
            if not (rho < 1):
                return arb('nan')
            # sup of F on [Y, Y+1]:  e^{a(Y+1)} P4 (2+Y)^4 * 2 (t+1) phi(t) with t >= Tb, (t+1)phi(t) decreasing for t >= 1
            Fsup = (aa * (Yb + 1)).exp() * P4 * (2 + Yb) ** 4 * 2 * (Tb + 1) * phi(Tb)
            ytail += mass * Exi * Fsup / (1 - rho)
        # z-tail |z| > L : envelope of the nu-form integrand (same quantity), as in Igen
        X = self.xi0 + self.xi1 * L
        e0, e1, e2 = self._envelope(2, X, lam_abs, w_lo, s_abs)
        PL = e0 + e1 * (X + arb('1.6')) + e2 * (1 + X * X + arb('0.8') * X)
        ztail = PL / (L * L) * T2(L)
        if verbose:
            print('I2_y', lam, 'core', tot, 'ytail', ytail, 'ztail', ztail, 'evals', nevals,
                  round(time.time() - t0, 1), flush=True)
        return self.a * (tot + arb(0, ub(ytail) + ub(ztail)))

    def Ip01_y(self, lam, s, inner_tol=None, outer_tol=None, L=None, verbose=False):
        """(I_s(lam), I_s'(lam)) at a point lam > 0 via the boundary-layer variable nu = xi + w y
        (see I2_y):  I = alpha E_z E(xi) int_0^inf K dy,  K = w e^Q h(x),
        K' = e^Q [ (w' + w Q') h + w x' h' ].  Packed as one nested acb.integral of K + i K'
        (both real on the real axis).  Tails: y > Y via the bound of I2_y with
        |h|, |h'| <= 2 (t+1) phi(t); |z| > L via the nu-form envelopes (kinds 0, 1) of Igen."""
        lam = arb(lam)
        s = arb(s)
        assert lam > 0
        w = (1 - lam * lam).sqrt()
        k, rq, s1q, kk = self.k, self.rq, self.w1q, self.kk
        w1 = -lam / w
        f = (1 - lam) / (1 + lam)
        r = f.sqrt()
        f1 = -2 / (1 + lam) ** 2
        r1 = f1 / (2 * r)
        lam_lo = arb(lb(lam))
        itol = inner_tol or self.tol
        otol = outer_tol or self.tol
        lam_c, w_c, w1_c, r_c, r1_c = acb(lam), acb(w), acb(w1), acb(r), acb(r1)
        s_c, kk_c = acb(s), acb(kk)
        I_ = acb(0, 1)

        def make_outer(Y):
            Yc = acb(Y)

            def outer(z, an_out):
                xi = (k - rq * z) / s1q
                Exi = 1 / (SQ2PI * Ps(xi))
                c = Exi * s_c / kk_c

                def inner(t, an_in):
                    an = an_in or an_out
                    y = Yc * t
                    Q = -(xi * w_c * y) - (w_c * y) ** 2 / 2
                    Q1 = -(xi * w1_c * y) - w_c * w1_c * y * y
                    x = xi * r_c - lam_c * y - c
                    x1 = xi * r1_c - y
                    h0 = lnPsi(x, an)
                    h1 = -Ecal(x)
                    eQ = Q.exp()
                    return eQ * (w_c * h0 + I_ * ((w1_c + w_c * Q1) * h0 + w_c * x1 * h1))

                if ub(z.rad()) < 1e-12:
                    J = acb.integral(inner, 0, 1, rel_tol=itol, abs_tol=itol * 1e-3)
                else:
                    J = acb.integral(inner, 0, 1, rel_tol=1e-4, abs_tol=itol, eval_limit=self.ball_eval_limit)
                return J * Yc * Exi * phi(z)
            return outer

        L = L or self.L
        tot = acb(0)
        ytail = arb(0)
        lam_abs = arb(ub(abs(lam)))
        w_lo = arb(lb(w))
        s_abs = arb(ub(abs(s)))
        t0 = time.time()
        for i in range(2 * L):
            za = arb(-L + i)
            zb = za + 1
            xa = (k - rq * za) / s1q
            xb = (k - rq * zb) / s1q
            xiR = xa.union(xb)
            xlo, xhi = lb(xiR), ub(xiR)
            ExiR = Ecal(acb(arb(xlo))).real.union(Ecal(acb(arb(xhi))).real)
            cR = ExiR * s / kk
            a_ = max(0.0, -xlo) * ub(w)
            b_ = lb(cR - xiR * r)
            T = max(9.0, (2.0 * a_ + 2.0) / float(lam_lo) + 9.0)
            Y = max(1.0, math.ceil((T - b_) / float(lam_lo)))
            _t = time.time()
            tot += acb.integral(make_outer(Y), za, zb, rel_tol=otol, abs_tol=otol * 1e-3)
            if verbose:
                print('piece', i, round(time.time() - _t, 2), tot, flush=True)
            # tail y > Y (same majorant as in I2_y; P(y) <= Pall (1+y)^4 covers both K and K')
            X = arb(max(abs(xlo), abs(xhi)))
            aw1, ar1 = abs(w1), abs(r1)
            cQ1 = X * aw1 + w * aw1
            cx1 = X * ar1 + 1
            Pall = w + (aw1 + w * cQ1) + w * cx1
            Yb = arb(Y)
            Tb = lam_lo * Yb + arb(b_)
            if not (Tb >= T - 1e-9):
                return None
            aa = arb(a_)
            rho = (1 + 1 / (1 + Yb)) ** 5 * (1 + lam_abs / (Tb + 1)) * (aa - lam_lo * lam_lo / 2 - lam_lo * Tb).exp()
            if not (rho < 1):
                return None
            Fsup = (aa * (Yb + 1)).exp() * Pall * (2 + Yb) ** 4 * 2 * (Tb + 1) * phi(Tb)
            ytail += (Psi(za) - Psi(zb)) * ExiR * Fsup / (1 - rho)
        X = self.xi0 + self.xi1 * L
        res = []
        for j, kd in enumerate([0, 1]):
            e0, e1, e2 = self._envelope(kd, X, lam_abs, w_lo, s_abs)
            PL = e0 + e1 * (X + arb('1.6')) + e2 * (1 + X * X + arb('0.8') * X)
            ztail = PL / (L * L) * T2(L)
            part = tot.real if j == 0 else tot.imag
            res.append(self.a * (part + arb(0, ub(ytail) + ub(ztail))))
        if verbose:
            print('Ip01_y', lam, res, 'ytail', ytail, round(time.time() - t0, 1), flush=True)
        return tuple(res)

    # ------------------------------------------------------------------ assembled pieces
    def PB_parts(self, lam, s, kinds=(0, 1), **kw):
        lam = arb(lam)
        s = arb(s)
        pq = self.psi * (1 - self.q)
        f = (1 - lam) / (1 + lam)
        sf = f.sqrt()
        fp = -2 / (1 + lam) ** 2
        fpp = 4 / (1 + lam) ** 3
        out = {}
        if 0 in kinds and 1 in kinds and kw.pop('ycoord', False):
            I0, I1 = self.Ip01_y(lam, s, **kw)
        elif 0 in kinds and 1 in kinds:
            I0, I1 = self.Igen('p01', lam, s, **kw)
        else:
            I0 = self.Igen(0, lam, s, **kw) if 0 in kinds else None
            I1 = self.Igen(1, lam, s, **kw) if 1 in kinds else None
        if 0 in kinds:
            out[0] = -self.Pstar() + pq * f / 2 + s * s / 2 - self.kk * s * sf + I0
        if 1 in kinds:
            out[1] = pq * fp / 2 - self.kk * s * fp / (2 * sf) + I1
        if 2 in kinds:
            I2 = self.Igen(2, lam, s, **kw)
            out[2] = pq * fpp / 2 - self.kk * s * (fpp / (2 * sf) - fp * fp / (4 * sf ** 3)) + I2
        return out

    def PB2_range(self, lam, s, hz=0.05, L=6):
        """Enclosure of PB''(lam) for every lam in the real ball lam (crude, for the second-order
        remainder): range enclosures I2_real (lam < 0.6) or I2_y (lam >= 0.6)."""
        lam = arb(lam)
        s = arb(s)
        pq = self.psi * (1 - self.q)
        f = (1 - lam) / (1 + lam)
        sf = f.sqrt()
        fp = -2 / (1 + lam) ** 2
        fpp = 4 / (1 + lam) ** 3
        I2 = self.I2_y(lam, s, hz=hz, L=L) if lb(lam) >= 0.6 else self.I2_real(lam, s, hz=hz, L=L)
        return pq * fpp / 2 - self.kk * s * (fpp / (2 * sf) - fp * fp / (4 * sf ** 3)) + I2
