"""Run one TD-MS parameter set through every animal and score it."""
from __future__ import annotations

import os
import pickle
import time

import numpy as np

import evaluate as E
import jeong as J
import sim as S
import sim_rich as R

GROUPS_MAIN = ["T12", "T3", "T4", "T6", "T5", "T7", "T8", "T9"]


class Context:
    def __init__(self, cache, ctx_file, fill_gaps=False, seq=True):
        self.cache = cache
        subjects = J.RR8 + (J.SEQ_CTRL + J.SEQ_INH if seq else [])
        self.animals = {s: S.Animal(cache, s, fill_gaps=fill_gaps) for s in subjects}
        self.sess = {s: a.sess for s, a in self.animals.items()}
        if os.path.exists(ctx_file):
            d = pickle.load(open(ctx_file, "rb"))
        else:
            used = set()

            def get(subj, name):
                used.add((subj, name))
                return self.animals[subj].raw[name]
            groups = GROUPS_MAIN + (["SEQ"] if seq else [])
            Md = E.measures_for(get, self.sess, groups)
            d = dict(data_stats=E.stats_of(Md), sds=E.noise_sds(Md),
                     refs={s: E.scale_ref(get, self.sess[s], s) for s in subjects},
                     inh_target={s: E.T.stat_T11({s: Md["SEQ_INH"][s]})["inhibition_ratio"][0]
                                 for s in J.SEQ_INH} if seq else {})
            d["used"] = used
            pickle.dump(d, open(ctx_file, "wb"))
        self.__dict__.update(d)


def simulate(ctx, p, subjects, inh=None):
    cls = p.get("cls", "MS")
    table = None if cls == "BS" else S.ms_table(p["m"], p["sigma"], p["d"], csc=p.get("csc", False))
    model = {}
    for s in subjects:
        I = (inh or {}).get(s, 0.0)
        if cls == "BS":
            out, _ = R.run_bs(ctx.animals[s], p, inh=I)
        elif cls == "SR":
            out, _ = R.run_sr(ctx.animals[s], p, table, inh=I)
        else:
            out, _ = ctx.animals[s].run(p, table, inh=I)
        for n, (st, delta) in out.items():
            if not np.all(np.isfinite(delta)):
                raise FloatingPointError("diverged")
        out = {n: v for n, v in out.items() if (s, n) in ctx.used}
        model[s] = ctx.animals[s].model_sessions(out)
    return model


def calibrate_inh(ctx, p, iters=3):
    """Per inhibited animal, I such that the model's inhibition ratio equals the data's."""
    import tests as T
    res = {}
    for s in J.SEQ_INH:
        target = ctx.inh_target[s]

        def ratio(I):
            m = simulate(ctx, p, [s], {s: I})[s]
            meas = T.measure_seq(lambda subj, name: m[name], ctx.sess[s], s)
            return meas["da_inh_cs2"] / meas["da_inh_rw"]
        I0, I1 = 0.0, 0.1
        r0, r1 = ratio(I0), ratio(I1)
        for _ in range(iters):
            if r1 == r0:
                break
            I2 = max(0.0, I1 + (target - r1) * (I1 - I0) / (r1 - r0))
            I0, r0, I1, r1 = I1, r1, I2, ratio(I2)
        res[s] = I1
    return res


def score(ctx, p, seq=True, nseeds=E.NSEEDS, inh=None):
    t0 = time.time()
    subjects = list(J.RR8)
    if seq:
        inh = inh if inh is not None else calibrate_inh(ctx, p)
        subjects += J.SEQ_CTRL + J.SEQ_INH
    model = simulate(ctx, p, subjects, inh)
    get = lambda subj, name: model[subj][name]
    groups = GROUPS_MAIN + (["SEQ"] if seq else [])
    Mm = E.measures_for(get, ctx.sess, groups)
    refs_m = {s: E.scale_ref(get, ctx.sess[s], s) for s in subjects}
    scales = {s: (ctx.refs[s] / refs_m[s]) if refs_m[s] > 0 else np.nan for s in subjects}
    st_seeds = [E.stats_of(E.add_noise(Mm, scales, ctx.sds, seed)) for seed in range(nseeds)]
    st = E.median_stats(st_seeds)
    P = E.passes(st, ctx.data_stats)
    mg = E.margins(st, ctx.data_stats)
    loss = float(sum(min(0.0, v) ** 2 if np.isfinite(v) else 100.0 for v in mg.values()))
    return dict(stats=st, passes=P, n_pass=int(sum(P.values())), margins=mg, loss=loss, scales=scales, inh=inh,
                seconds=time.time() - t0)
