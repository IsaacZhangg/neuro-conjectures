# AGENT-WRITTEN (Claude), NOT INDEPENDENTLY REVIEWED.
#
# Box-by-box certification of  S(lambda) < 0  on  [lam_min, -delta] U [delta, 1-eps]
# for the Ding--Sun second-moment exponent at margin kappa (see dsarb.py, README.md).
#
# For each box [A0, A1] in the overlap-profile parameter A (lambda = ell(A) increasing), with an
# exact interior point lam* and a fixed s:
#   S(lam) <= H(lam) + PB(lam; s)                                        (A(lam) <= B(lam, s))
#   H(lam) <= H(Ac) + H'(lam_c) (lam - lam_c),   lam_c = ell(Ac),  H'(lam_c) = -(1-q) ln(Ac)/2
#                                                                        (H concave in lam)
#   PB(lam) <= PB(lam*) + PB'(lam*) d + (1/2) max(0, sup_Lambda PB'') d^2,    d = lam - lam*
# so with g = H'(lam_c) + PB'(lam*),  C = H(Ac) + H'(lam_c)(lam* - lam_c) + PB(lam*):
#   sup_box S <= C + max_{d in {d0, d1}} [ g d + M d^2 / 2 ]     (convex in d -> endpoint max)
# where d0 = lb(ell(A0)) - lam*, d1 = ub(ell(A1)) - lam*, Lambda = [lb ell(A0), ub ell(A1)].
# All quantities are Arb enclosures valid for every (alpha, q, psi) in the certified parameter box.
import sys, os, json, time, math
import multiprocessing as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
FLOATDIR = os.environ.get('DS_FLOAT_DIR', '/tmp/claude-0/-home-user/fa4b6461-5216-5a3d-ae79-fb39767abede/scratchpad/dingsun')
sys.path.insert(0, FLOATDIR)

PREC = 96
TOL = 1e-11
MW = float(os.environ.get('DS_MW', '0.004'))
MWMAX = float(os.environ.get('DS_MWMAX', '0.02'))
LM = int(os.environ.get('DS_LM', '6'))   # z-truncation for the (crude) sup PB'' enclosure
HZ = float(os.environ.get('DS_HZ', '0.05'))  # z-cell width of the range enclosure of PB''

# ---------------------------------------------------------------- float design (NOT rigorous)
def float_model(kappa, alpha):
    import ds_S as d
    return d.DS(float(kappa), alpha=float(alpha))


def f_PB(ds, lam, s):
    return ds.Pfun(lam) + ds.Bfun(lam, s)


def f_est(ds, lam):
    """float: s*(lam), Sfix, Sfix', PB'' (s fixed)."""
    h, s = ds.Hfun(lam), ds.Afun(lam, True)[1]
    S0 = h + f_PB(ds, lam, s)
    e = 1e-5 * min(1.0, 1 - lam, abs(lam) + 1e-3) + 1e-9
    Sp = (ds.Hfun(lam + e) + f_PB(ds, lam + e, s) - ds.Hfun(lam - e) - f_PB(ds, lam - e, s)) / (2 * e)
    e2 = 1e-3 * min(1.0, 1 - lam)
    P2 = (f_PB(ds, lam + e2, s) - 2 * f_PB(ds, lam, s) + f_PB(ds, lam - e2, s)) / e2 ** 2
    return s, S0, Sp, P2


def design(ds, lo, hi, frac=0.45, rmax=0.04, first_at_lammin=False):
    """greedy float design of lambda-boxes covering [lo, hi]; returns list of (la, lstar, lb, s)."""
    boxes = []
    la = lo
    while la < hi:
        r = min(rmax, 0.25 * (1 - la), 0.3 * max(abs(la), 1e-3))
        for _ in range(4):
            ls = min(la + r, hi)
            s, S0, Sp, P2 = f_est(ds, ls)
            Mx = max(P2, 0) + 0.5 * abs(P2) + 1.0   # allowance for the crude PB'' enclosure
            # want |Sp| r + Mx r^2/2 <= frac |S0|
            a, b, c = Mx / 2, abs(Sp), -frac * abs(S0)
            rn = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
            rn = min(rn, rmax, 0.6 * (1 - la))
            if abs(rn - r) < 0.05 * r:
                r = rn
                break
            r = rn
        ls = la + r
        lb_ = la + 2 * r
        if lb_ >= hi:
            lb_ = hi
            ls = 0.5 * (la + hi)
        s = f_est(ds, ls)[0]
        boxes.append([la, ls, lb_, round(s, 7)])
        la = lb_
    return boxes


# ---------------------------------------------------------------- rigorous part
_M = None


def _init(kappa, pfile):
    global _M
    from flint import ctx
    ctx.prec = PREC
    import dsarb
    P = json.load(open(pfile))
    from flint import arb
    al = dsarb.box(arb(P['alpha_lb']), arb(P['alpha_ub']))
    _M = dsarb.Model(kappa, al, arb(P['q0']), arb(P['psi0']), tol=TOL)
    _M.Hstar(); _M.Pstar()


def certify_box(spec):
    import dsarb
    from dsarb import arb, ub, lb, box
    M = _M
    t0 = time.time()
    A0, A1, Ac = spec['A0'], spec['A1'], spec['Ac']
    ls = arb(spec['lstar'])
    s = arb(spec['s'])
    l0 = M.lam_min() if float(A0) == 0.0 else M.ell(arb(A0))
    l1 = M.ell(arb(A1))
    lc = M.ell(arb(Ac))
    Hc = M.Hfun(arb(Ac))
    Hpc = M.Hprime(arb(Ac))
    t1 = time.time()
    PB = M.PB_parts(ls, s, kinds=(0, 1), ycoord=(float(spec['lstar']) >= 0.6))
    Lam = box(lb(l0), ub(l1))
    t2 = time.time()
    # sup of PB'' over Lambda: split into sub-balls of width <= MW (wide lambda-balls are slow)
    lo_, hi_ = lb(l0), ub(l1)
    # sub-balls of width <= min(MWMAX, 0.1 (1 - lambda)) (the integrand scale grows like (1-lam)^{-3/2})
    cuts = [lo_]
    while cuts[-1] < hi_:
        x = cuts[-1]
        step = max(min(MWMAX, 0.1 * (1.0 - min(hi_, x + MWMAX))), 2e-4)
        cuts.append(min(hi_, x + step))
    nsub = len(cuts) - 1
    PB2 = None
    for j in range(nsub):
        sub = box(cuts[j], cuts[j + 1])
        v = M.PB2_range(sub, s, hz=HZ, L=LM)
        PB2 = v if PB2 is None else (PB2 if ub(PB2) >= ub(v) else v)
    Mx = arb(max(ub(PB2), 0.0))
    t3 = time.time()
    g = Hpc + PB[1]
    C = Hc + Hpc * (ls - lc) + PB[0]
    d0 = arb(lb(l0)) - ls
    d1 = arb(ub(l1)) - ls
    f0 = g * d0 + Mx * d0 * d0 / 2
    f1 = g * d1 + Mx * d1 * d1 / 2
    bound = ub(C) + max(ub(f0), ub(f1))
    out = dict(spec)
    out.update(l0=str(l0), l1=str(l1), lc=str(lc), H_Ac=str(Hc), Hprime=str(Hpc),
               PB=str(PB[0]), PBprime=str(PB[1]), PB2_sup_enclosure=str(PB2),
               Slin_const=str(C), slope=str(g), bound=bound, ok=bool(bound < 0),
               l0_lb=lb(l0), l1_ub=ub(l1), l0_ub=ub(l0), l1_lb=lb(l1),
               seconds=round(time.time() - t0, 1), t_1d=round(t1 - t0, 1), t_PB01=round(t2 - t1, 1),
               t_PB2=round(t3 - t2, 1), n_sub=nsub)
    return out


def make_specs(ds, boxes, A_left=None):
    specs = []
    prevA = A_left
    for i, (la, lsv, lbv, s) in enumerate(boxes):
        A0 = prevA if prevA is not None else repr(float(ds.A_of_lam(la)))
        A1 = repr(float(ds.A_of_lam(lbv)))
        Ac = repr(float(ds.A_of_lam(lsv)))
        specs.append(dict(A0=A0, A1=A1, Ac=Ac, lstar=repr(float(lsv)), s=repr(float(s)), lam_float=[float(la), float(lsv), float(lbv)]))
        prevA = A1
    return specs


def bisect_spec(ds, sp):
    la, lsv, lbv = sp['lam_float']
    lm = 0.5 * (la + lbv)
    Am = repr(float(ds.A_of_lam(lm)))
    out = []
    for (a, b, A0, A1) in [(la, lm, sp['A0'], Am), (lm, lbv, Am, sp['A1'])]:
        c = 0.5 * (a + b)
        s = round(ds.Afun(c, True)[1], 7)
        out.append(dict(A0=A0, A1=A1, Ac=repr(float(ds.A_of_lam(c))), lstar=repr(float(c)), s=repr(float(s)), lam_float=[a, c, b],
                        parent=sp.get('id')))
    return out


def run(kappa, pfile, delta, eps, outdir, nproc=2, max_rounds=6, tag=''):
    os.makedirs(outdir, exist_ok=True)
    P = json.load(open(pfile))
    from flint import arb
    amid = float((arb(P['alpha_lb']) + arb(P['alpha_ub'])) / 2)
    ds = float_model(kappa, amid)
    lmin = ds.lam_min
    t_start = time.time()
    neg = design(ds, lmin, -delta + 1e-7)
    pos = design(ds, delta - 1e-7, 1 - eps + 1e-7)
    specs = make_specs(ds, neg, A_left='0.0') + make_specs(ds, pos)
    for i, sp in enumerate(specs):
        sp['id'] = str(i)
    log = open(os.path.join(outdir, f'log_kappa{kappa}{tag}.txt'), 'a')
    log.write(f'# start {time.ctime()} kappa={kappa} delta={delta} eps={eps} nboxes={len(specs)}\n'); log.flush()
    done = []
    pending = specs
    rnd = 0
    with mp.Pool(nproc, initializer=_init, initargs=(kappa, pfile)) as pool:
        while pending and rnd < max_rounds:
            nxt = []
            for r in pool.imap_unordered(certify_box, pending):
                r['round'] = rnd
                log.write(json.dumps({k: r[k] for k in ('id', 'lam_float', 'bound', 'ok', 'seconds')}) + '\n'); log.flush()
                if r['ok']:
                    done.append(r)
                else:
                    r['failed'] = True
                    done.append(r)
                    kids = bisect_spec(ds, r)
                    for j, k in enumerate(kids):
                        k['id'] = r['id'] + '.' + str(j)
                    nxt += kids
            pending = nxt
            rnd += 1
    # ---------------- assemble certificate: the leaves that passed, and a coverage check
    leaves = [r for r in done if r['ok']]
    unresolved = pending
    failed_leaves = [r for r in done if not r['ok']]
    leaves.sort(key=lambda r: r['lam_float'][0])
    cert = dict(kappa=kappa, params=P, delta=delta, eps=eps, prec=PREC, tol=TOL,
                n_boxes_passed=len(leaves), n_unresolved=len(unresolved),
                wall_seconds=round(time.time() - t_start, 1))
    # coverage: passed leaves must tile each side contiguously in A
    def chain(rs, A_first):
        Acur = A_first
        for r in rs:
            if r['A0'] != Acur:
                return False, Acur
            Acur = r['A1']
        return True, Acur
    negl = [r for r in leaves if r['lam_float'][2] <= 0]
    posl = [r for r in leaves if r['lam_float'][0] >= 0]
    okn, An = chain(negl, '0.0')
    okp, Ap = (chain(posl, posl[0]['A0']) if posl else (False, None))
    cert['coverage'] = dict(
        negative_chain_contiguous=okn, negative_last_A=An,
        negative_right_end_lb=(negl[-1]['l1_lb'] if negl else None),
        negative_reaches_minus_delta=bool(negl and negl[-1]['l1_lb'] >= -delta),
        positive_chain_contiguous=okp,
        positive_left_end_ub=(posl[0]['l0_ub'] if posl else None),
        positive_starts_below_delta=bool(posl and posl[0]['l0_ub'] <= delta),
        positive_right_end_lb=(posl[-1]['l1_lb'] if posl else None),
        positive_reaches_1_minus_eps=bool(posl and posl[-1]['l1_lb'] >= 1 - eps))
    cert['certified'] = bool(okn and okp and not unresolved and cert['coverage']['negative_reaches_minus_delta']
                             and cert['coverage']['positive_starts_below_delta']
                             and cert['coverage']['positive_reaches_1_minus_eps'])
    worst = max(leaves, key=lambda r: r['bound']) if leaves else None
    cert['worst_bound'] = worst['bound'] if worst else None
    cert['worst_box_lambda'] = worst['lam_float'] if worst else None
    cert['boxes'] = leaves
    cert['failed_then_bisected'] = failed_leaves
    cert['unresolved'] = unresolved
    fn = os.path.join(outdir, f'certificate_kappa{kappa}{tag}.json')
    json.dump(cert, open(fn, 'w'), indent=1)
    log.write(f'# end {time.ctime()} certified={cert["certified"]} worst={cert["worst_bound"]} at {cert["worst_box_lambda"]}\n')
    log.close()
    return cert


if __name__ == '__main__':
    kappa = sys.argv[1]
    pfile = sys.argv[2]
    delta = float(sys.argv[3]); eps = float(sys.argv[4])
    outdir = sys.argv[5] if len(sys.argv) > 5 else os.path.join(HERE, 'out')
    tag = sys.argv[6] if len(sys.argv) > 6 else ''
    c = run(kappa, pfile, delta, eps, outdir, tag=tag)
    print(json.dumps({k: v for k, v in c.items() if k not in ('boxes', 'failed_then_bisected', 'unresolved', 'params')}, indent=1))
