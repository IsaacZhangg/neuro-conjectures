"""Read only the event log and dF/F of the sequential-conditioning (stGtACR) sessions
from DANDI 000351 by HTTP range requests, without downloading raw photometry.

Usage: python fetch_opto.py ASSETS_YAML CACHE_DIR
ASSETS_YAML is the dandiset version's asset list (path, contentUrl, contentSize).
"""
import json, os, sys
import fsspec, h5py, numpy as np, yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract import extract_h5  # noqa: E402


def main(assets, cache):
    a = yaml.safe_load(open(assets))
    fs = fsspec.filesystem("https", block_size=2 ** 21)
    man = {}
    for x in a:
        p = x["path"]
        if "stGtACR" not in p:
            continue
        subj = p.split("/")[0][4:]
        ses = p.split("_ses-")[1][:-4]
        dst = os.path.join(cache, subj, ses + ".npz")
        man[p] = dict(size=x["contentSize"], url=[u for u in x["contentUrl"] if "s3" in u][0])
        if os.path.exists(dst):
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        with fs.open(man[p]["url"], "rb") as f, h5py.File(f, "r") as h:
            np.savez_compressed(dst, **extract_h5(h))
        print(p, flush=True)
    json.dump(man, open(os.path.join(cache, "manifest_opto.json"), "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
