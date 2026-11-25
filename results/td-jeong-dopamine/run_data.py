"""Compute Tests 1 to 11 on the real dopamine data. Usage: python run_data.py CACHE OUT.json"""
import json, sys
import numpy as np
import jeong as J
import tests as T


def collect(cache):
    sess = {s: J.list_sessions(cache, s) for s in J.RR8 + T.__dict__.get("EXTRA", [])}
    cachemap = {}

    def get(subj, name):
        k = (subj, name)
        if k not in cachemap:
            cachemap.clear()
            cachemap[k] = J.load(cache, subj, name)
        return cachemap[k]
    return sess, get


def measures(cache, subjects_main=J.MAIN7, rr=J.RR8, seq=True):
    sess = {s: J.list_sessions(cache, s) for s in rr + (J.SEQ_CTRL + J.SEQ_INH if seq else [])}
    get = lambda subj, name: J.load(cache, subj, name)
    M = {}
    M["T12"] = {s: T.measure_randomrewards(get, sess[s], s) for s in rr}
    for iC, key in ((0, "T3"), (1, "T4"), (2, "T6")):
        M[key] = {s: T.measure_cumsum(get, sess[s], s, iC) for s in subjects_main}
    M["T5"] = {s: T.measure_T5(get, sess[s], s) for s in subjects_main}
    M["T7"] = {s: T.measure_T7(get, sess[s], s) for s in subjects_main}
    M["T8"] = {s: T.measure_T8(get, sess[s], s) for s in subjects_main}
    M["T9"] = {s: T.measure_T9(get, sess[s], s) for s in subjects_main}
    if seq:
        M["SEQ_CTRL"] = {s: T.measure_seq(get, sess[s], s) for s in J.SEQ_CTRL}
        M["SEQ_INH"] = {s: T.measure_seq(get, sess[s], s) for s in J.SEQ_INH}
    return M


def all_stats(M):
    out = dict(T1=T.stat_T1(M["T12"]), T2=T.stat_T2(M["T12"]), T3=T.stat_T3(M["T3"]), T4=T.stat_T4(M["T4"]),
               T5=T.stat_T5(M["T5"]), T6=T.stat_T6(M["T6"]), T7=T.stat_T7(M["T7"]), T8=T.stat_T8(M["T8"]),
               T9=T.stat_T9(M["T9"]))
    if "SEQ_CTRL" in M:
        out["T10"] = T.stat_T10(M["SEQ_CTRL"])
        out["T11"] = T.stat_T11(M["SEQ_INH"])
    return out


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


if __name__ == "__main__":
    seq = "--noseq" not in sys.argv
    M = measures(sys.argv[1], seq=seq)
    st = all_stats(M)
    json.dump(clean(st), open(sys.argv[2], "w"), indent=1)
    for k, v in st.items():
        print(k, {kk: (round(vv["t"], 2), round(vv["p"], 4)) if isinstance(vv, dict) and "t" in vv else None
                  for kk, vv in v.items()} if "t" not in v else (round(v["t"], 2), round(v["p"], 4), round(v["mean"], 3)))
