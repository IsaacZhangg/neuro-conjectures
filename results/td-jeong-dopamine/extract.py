"""Extract event logs and dF/F from local DANDI 000351 NWB files into .npz cache.

Usage: python extract.py RAW_DIR CACHE_DIR
Raw files are named <subject>__<subject>_ses-<session>.nwb (path with / -> __).
Writes CACHE_DIR/<subject>/<session>.npz and CACHE_DIR/manifest.json (size, sha256).
"""
import hashlib, json, os, re, sys
import h5py
import numpy as np


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def extract_h5(h):
    e = h["acquisition/eventlog"]
    out = dict(
        eventindex=e["eventindex"][:].astype(np.int32),
        eventtime=e["eventtime"][:].astype(np.float64),
        nonsolenoidflag=e["nonsolenoidflag"][:].astype(np.int8),
        desc=np.array([x.decode() for x in h["acquisition/eventidx_table/event_description"][:]]),
        code=h["acquisition/eventidx_table/event_index"][:].astype(np.int32),
        session_description=np.array(h["session_description"][()].decode()),
    )
    p = h["processing/photometry/dff"]
    out["dff"] = p["data"][:].astype(np.float32)
    out["dff_t"] = p["timestamps"][:].astype(np.float64)
    return out


def main(raw, cache):
    man = {}
    for fn in sorted(os.listdir(raw)):
        m = re.match(r"sub-(.+?)__sub-.+?_ses-(.+)\.nwb$", fn)
        if not m:
            continue
        subj, ses = m.groups()
        os.makedirs(os.path.join(cache, subj), exist_ok=True)
        dst = os.path.join(cache, subj, ses + ".npz")
        path = os.path.join(raw, fn)
        if not os.path.exists(dst):
            with h5py.File(path, "r") as h:
                np.savez_compressed(dst, **extract_h5(h))
        man[fn] = dict(size=os.path.getsize(path), sha256=sha256(path))
    json.dump(man, open(os.path.join(cache, "manifest.json"), "w"), indent=1)
    print(len(man), "files")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
