"""Reproducible check: Rall's 3/2 power rule at dendritic bifurcations.

Claim under test (Rall 1959; Scholarpedia "Rall model", 2009): at a dendritic
bifurcation, d0^(3/2) = d1^(3/2) + d2^(3/2). We estimate the branch exponent
eta in d0^eta = d1^eta + d2^eta per cell from NeuroMorpho.org reconstructions
and apply the decision rule frozen in CLAIM.md and METHODS.md.

Run:  python check.py            full run: downloads NeuroMorpho.org metadata and
                                 CNG SWC files into ./data/ (not committed),
                                 writes check-output.json and manifest.json
      python check.py --offline  full analysis from ./data/ only (no network)
      python check.py --bluebrain  secondary test (METHODS.md, Amendment A1) on
                                 Blue Brain Open Data; writes
                                 check-output-bluebrain.json
      python check.py --quick    estimator validation on synthetic trees only
                                 (CI smoke test; does not set the status)
Requires Python 3.10+ and NumPy. Deterministic given the seed and the data.
"""
from __future__ import annotations

import argparse
import hashlib
import http.client
import json
import platform
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

SEED = 20261125
RALL_ETA = 1.5
MARGIN = (1.35, 1.65)
ETA_GRID = np.round(np.arange(0.5, 6.0 + 1e-9, 0.005), 3)
N_BOOT = 2000
WINDOW_UM = 5.0
MIN_GROUP_ELIGIBLE = 10
MAX_PER_GROUP = 20
MIN_GROUP_QUALIFYING = 10
MIN_DEND_SAMPLES = 50
MIN_DISTINCT_DIAMS = 10
MAX_MODAL_FRACTION = 0.5
MIN_BIFURCATIONS = 10
MIN_GROUPS_FOR_VERDICT = 5
DENDRITE_TYPES = (3, 4)

API = "https://neuromorpho.org/api/neuron"
SWC_URL = "https://neuromorpho.org/dableFiles/{archive}/CNG%20version/{name}.CNG.swc"
USER_AGENT = "neuro-conjectures-check/1.0"
HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

# ---------------------------------------------------------------- SWC parsing


def parse_swc(text: str) -> dict:
    """Return arrays type, xyz, diameter, parent index (-1 for roots)."""
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        rows.append((int(float(parts[0])), int(float(parts[1])), float(parts[2]),
                     float(parts[3]), float(parts[4]), float(parts[5]), int(float(parts[6]))))
    if not rows:
        raise ValueError("no SWC rows")
    ids = [r[0] for r in rows]
    index = {i: k for k, i in enumerate(ids)}
    parent = np.array([index.get(r[6], -1) for r in rows])
    return {
        "type": np.array([r[1] for r in rows]),
        "xyz": np.array([[r[2], r[3], r[4]] for r in rows]),
        "diam": 2.0 * np.array([r[5] for r in rows]),
        "parent": parent,
    }


class Tree:
    def __init__(self, swc: dict):
        self.type = swc["type"]
        self.xyz = swc["xyz"]
        self.diam = swc["diam"]
        self.parent = swc["parent"]
        n = len(self.type)
        self.children: list[list[int]] = [[] for _ in range(n)]
        for k, p in enumerate(self.parent):
            if p >= 0:
                self.children[p].append(k)
        self.seglen = np.zeros(n)
        has_p = self.parent >= 0
        self.seglen[has_p] = np.linalg.norm(self.xyz[has_p] - self.xyz[self.parent[has_p]], axis=1)
        self.is_dend = np.isin(self.type, DENDRITE_TYPES)

    def dend_children(self, k: int) -> list[int]:
        return [c for c in self.children[k] if self.is_dend[c]]

    def is_branch(self, k: int) -> bool:
        return len(self.children[k]) >= 2

    def bifurcations(self) -> list[int]:
        out = []
        for k in range(len(self.type)):
            if not self.is_dend[k]:
                continue
            p = self.parent[k]
            if p < 0 or not self.is_dend[p]:
                continue
            if len(self.children[k]) == 2 and len(self.dend_children(k)) == 2:
                out.append(k)
        return out

    def parent_side(self, k: int, window: float | None) -> list[int]:
        """Node k and ancestors within `window` um, stopping before a branch point
        or a non-dendritic node. window=None means the whole section."""
        nodes = [k]
        dist = 0.0
        cur = k
        while True:
            p = self.parent[cur]
            if p < 0 or not self.is_dend[p] or self.is_branch(p):
                break
            dist += self.seglen[cur]
            if window is not None and dist > window:
                break
            nodes.append(p)
            cur = p
        return nodes

    def daughter_side(self, c: int, window: float | None) -> list[int]:
        """Samples from child c distally, within `window` um of the bifurcation,
        up to the next branch point or terminal. Always includes c."""
        nodes = [c]
        dist = self.seglen[c]
        cur = c
        while len(self.children[cur]) == 1:
            nxt = self.children[cur][0]
            if not self.is_dend[nxt]:
                break
            dist += self.seglen[nxt]
            if window is not None and dist > window:
                break
            nodes.append(nxt)
            cur = nxt
        return nodes

    def order(self, k: int) -> int:
        """Centrifugal order: 1 + number of dendritic bifurcations proximal to k."""
        n = 1
        p = self.parent[k]
        while p >= 0 and self.is_dend[p]:
            if len(self.dend_children(p)) >= 2:
                n += 1
            p = self.parent[p]
        return n


def branch_triples(tree: Tree, window: float | None) -> tuple[np.ndarray, np.ndarray]:
    """Rows (d0, d1, d2) for valid bifurcations, and their centrifugal orders."""
    rows, orders = [], []
    for k in tree.bifurcations():
        c1, c2 = tree.dend_children(k)
        d0 = tree.diam[tree.parent_side(k, window)].mean()
        d1 = tree.diam[tree.daughter_side(c1, window)].mean()
        d2 = tree.diam[tree.daughter_side(c2, window)].mean()
        if d0 > 0 and d1 > 0 and d2 > 0:
            rows.append((d0, d1, d2))
            orders.append(tree.order(k))
    return np.array(rows).reshape(-1, 3), np.array(orders, dtype=int)


# ---------------------------------------------------------------- estimator


def fit_eta(triples: np.ndarray) -> float:
    """eta in ETA_GRID minimising sum_b [(1/eta) ln(d1^eta + d2^eta) - ln d0]^2."""
    logd = np.log(triples)
    eta = ETA_GRID[:, None]
    pred = np.logaddexp(eta * logd[None, :, 1], eta * logd[None, :, 2]) / eta
    loss = ((pred - logd[None, :, 0]) ** 2).sum(axis=1)
    return float(ETA_GRID[int(np.argmin(loss))])


def within_10pct(triples: np.ndarray) -> float:
    lhs = triples[:, 0] ** 1.5
    rhs = triples[:, 1] ** 1.5 + triples[:, 2] ** 1.5
    return float(np.mean(np.abs(rhs / lhs - 1.0) <= 0.10))


def boot_median_ci(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    idx = rng.integers(0, len(values), size=(N_BOOT, len(values)))
    meds = np.median(values[idx], axis=1)
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def label(ci: tuple[float, float]) -> str:
    lo, hi = ci
    if lo >= MARGIN[0] and hi <= MARGIN[1]:
        return "consistent"
    if hi < MARGIN[0] or lo > MARGIN[1]:
        return "deviating"
    return "inconclusive"


# ---------------------------------------------------------------- synthetic


def synthetic_swc(eta: float, rng: np.random.Generator) -> str:
    """Binary dendritic tree obeying d0^eta = d1^eta + d2^eta, then 10% lognormal
    noise on every sample and 0.1 um quantisation."""
    lines = ["1 1 0 0 0 5 -1"]
    nid = 1

    def grow(parent_id: int, start: np.ndarray, direction: np.ndarray, d: float, depth: int):
        nonlocal nid
        n_samples = int(rng.integers(8, 25))
        pos = start.copy()
        last = parent_id
        for _ in range(n_samples):
            pos = pos + direction
            nid += 1
            dn = d * np.exp(rng.normal(0.0, 0.10))
            dn = max(0.1, round(dn / 0.1) * 0.1)
            lines.append(f"{nid} 3 {pos[0]:.3f} {pos[1]:.3f} {pos[2]:.3f} {dn / 2:.4f} {last}")
            last = nid
        if depth >= 7 or d < 0.5:
            return
        f = rng.uniform(0.2, 0.8)
        for share in (f, 1.0 - f):
            dd = d * share ** (1.0 / eta)
            nd = direction + rng.normal(0.0, 0.5, 3)
            grow(last, pos, nd / np.linalg.norm(nd), dd, depth + 1)

    for _ in range(int(rng.integers(3, 6))):
        v = rng.normal(size=3)
        grow(1, np.zeros(3), v / np.linalg.norm(v), float(rng.uniform(1.5, 3.5)), 0)
    return "\n".join(lines) + "\n"


def validate_estimator(n_cells: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    out = {}
    ok = True
    for true_eta in (1.5, 2.0):
        etas = []
        for _ in range(n_cells):
            tree = Tree(parse_swc(synthetic_swc(true_eta, rng)))
            triples, _ = branch_triples(tree, WINDOW_UM)
            etas.append(fit_eta(triples))
        med = float(np.median(etas))
        passed = abs(med - true_eta) <= 0.10
        ok &= passed
        out[str(true_eta)] = {"median_eta_hat": round(med, 3), "bias": round(med - true_eta, 3),
                              "cells": n_cells, "pass": bool(passed)}
    out["pass"] = bool(ok)
    return out


# ---------------------------------------------------------------- NeuroMorpho


def fetch(url: str, tries: int = 6, max_wait: float = 120.0) -> bytes:
    """GET url, retrying server errors (5xx, 429), dropped connections and
    timeouts with exponential backoff (2, 4, 8, 16, 32 s, or the server's
    Retry-After). Other 4xx errors, such as 404, are raised at once."""
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    err: Exception | None = None
    for attempt in range(tries):
        wait = min(2.0 ** (attempt + 1), max_wait)
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code < 500 and e.code != 429:
                raise
            err = e
            retry_after = e.headers.get("Retry-After") if e.headers else None
            if retry_after and retry_after.isdigit():
                wait = min(float(retry_after), max_wait)
        except (urllib.error.URLError, http.client.HTTPException, TimeoutError,
                ConnectionError) as e:
            err = e
        if attempt + 1 < tries:
            print(f"fetch failed ({err}); retry {attempt + 1}/{tries - 1} in {wait:.0f} s",
                  file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"could not fetch {url} after {tries} tries: {err}")


def load_metadata(offline: bool) -> tuple[list[dict], dict]:
    path = DATA / "metadata.json"
    if path.exists():
        cached = json.loads(path.read_text())
        return cached["neurons"], cached["info"]
    if offline:
        raise SystemExit("no cached metadata in ./data; run without --offline first")
    # Each page is saved as it arrives, so a run that stops part way (for
    # example on a server error) resumes from the first missing page.
    keep = ["neuron_id", "neuron_name", "archive", "species", "cell_type", "brain_region",
            "attributes", "domain", "experiment_condition", "reference_pmid", "reference_doi"]
    pages_dir = DATA / "metadata-pages"
    pages_dir.mkdir(parents=True, exist_ok=True)
    page, total_pages, total = 0, None, None
    while total_pages is None or page < total_pages:
        page_file = pages_dir / f"page-{page:05d}.json"
        if page_file.exists():
            saved = json.loads(page_file.read_text())
        else:
            body = json.loads(fetch(f"{API}?page={page}&size=500&sort=neuron_id,asc"))
            saved = {"total_pages": body["page"]["totalPages"],
                     "total_elements": body["page"]["totalElements"],
                     "neurons": [{k: n.get(k) for k in keep}
                                 for n in body.get("_embedded", {}).get("neuronResources", [])]}
            tmp = page_file.with_suffix(".tmp")
            tmp.write_text(json.dumps(saved))
            tmp.replace(page_file)
            time.sleep(0.2)
        total_pages, total = saved["total_pages"], saved["total_elements"]
        page += 1
        print(f"metadata page {page}/{total_pages}", file=sys.stderr)
    neurons, seen = [], set()
    for page_file in sorted(pages_dir.glob("page-*.json"))[:total_pages]:
        for n in json.loads(page_file.read_text())["neurons"]:
            if n["neuron_id"] not in seen:
                seen.add(n["neuron_id"])
                neurons.append(n)
    info = {"download_date_utc": time.strftime("%Y-%m-%d", time.gmtime()),
            "api_total_elements": total, "neurons_listed": len(neurons)}
    DATA.mkdir(exist_ok=True)
    path.write_text(json.dumps({"info": info, "neurons": neurons}))
    return neurons, info


def is_control(cond) -> bool:
    """True only for control animals. The API gives experiment_condition as a
    list, such as ['Control']; a plain string is accepted too."""
    if isinstance(cond, str):
        return cond == "Control"
    return isinstance(cond, list) and len(cond) > 0 and all(c == "Control" for c in cond)


def eligible_by_metadata(n: dict) -> bool:
    attrs = n.get("attributes") or ""
    domain = n.get("domain") or ""
    cells = n.get("cell_type") or []
    return (attrs.startswith("Diameter")
            and "Dendrites" in domain and "No Dendrites" not in domain
            and is_control(n.get("experiment_condition"))
            and not any("glia" in str(c).lower() for c in cells))


def group_key(n: dict) -> str:
    return " | ".join([str(n.get("archive")), str(n.get("species")),
                       ", ".join(str(c) for c in (n.get("cell_type") or []))])


def load_swc(n: dict, offline: bool) -> tuple[str | None, str | None]:
    path = DATA / "swc" / f"{n['neuron_id']}.swc"
    if path.exists():
        return path.read_text(errors="replace"), None
    if offline:
        return None, "not cached"
    url = SWC_URL.format(archive=urllib.parse.quote(str(n["archive"]).lower()),
                         name=urllib.parse.quote(str(n["neuron_name"])))
    try:
        raw = fetch(url)
    except Exception as e:  # noqa: BLE001 - recorded as an exclusion reason
        return None, f"download failed: {type(e).__name__}"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)
    time.sleep(0.2)
    return raw.decode(errors="replace"), None


def file_filter(tree: Tree) -> str | None:
    d = tree.diam[tree.is_dend]
    if d.size < MIN_DEND_SAMPLES:
        return "fewer than 50 dendritic samples"
    vals, counts = np.unique(np.round(d, 4), return_counts=True)
    if vals.size < MIN_DISTINCT_DIAMS:
        return "fewer than 10 distinct diameters"
    if counts.max() / d.size > MAX_MODAL_FRACTION:
        return "modal diameter over 50% of samples"
    return None


# ---------------------------------------------------------------- analysis


VARIANT_WINDOWS = {"primary_5um": WINDOW_UM, "R1_section_mean": None, "R2_10um": 10.0}
ANALYSES = ["primary_5um", "R1_section_mean", "R2_10um", "R3_drop_identical"]


def cell_triples(text: str) -> tuple[dict, dict, str | None]:
    """Parse one SWC file and apply the file filters. Returns (triples, orders, exclusion reason)."""
    try:
        tree = Tree(parse_swc(text))
    except Exception:  # noqa: BLE001 - recorded as an exclusion reason
        return {}, {}, "parse failed"
    why = file_filter(tree)
    if why is None and len(tree.bifurcations()) < MIN_BIFURCATIONS:
        why = "fewer than 10 bifurcations"
    triples, orders = {}, {}
    if why is None:
        for name, w in VARIANT_WINDOWS.items():
            triples[name], orders[name] = branch_triples(tree, w)
        if len(triples["primary_5um"]) < MIN_BIFURCATIONS:
            why = "fewer than 10 valid bifurcations"
    if why is not None:
        return {}, {}, why
    triples["R3_drop_identical"] = triples["primary_5um"]
    orders["R3_drop_identical"] = orders["primary_5um"]
    return triples, orders, None


def analyse(cells: list[dict], variant: str, rng: np.random.Generator) -> dict:
    per_cell = []
    for c in cells:
        tr, orders = c["triples"][variant], c["orders"][variant]
        if variant == "R3_drop_identical":
            keep = ~((tr[:, 0] == tr[:, 1]) & (tr[:, 1] == tr[:, 2]))
            tr, orders = tr[keep], orders[keep]
        if len(tr) < MIN_BIFURCATIONS:
            continue
        per_cell.append({"group": c["group"], "archive": c["archive"], "eta": fit_eta(tr),
                         "within10": within_10pct(tr),
                         "eta_proximal": fit_eta(tr[orders <= 2]) if (orders <= 2).sum() >= 3 else None,
                         "eta_distal": fit_eta(tr[orders >= 3]) if (orders >= 3).sum() >= 3 else None})
    by_group = defaultdict(list)
    for c in per_cell:
        by_group[c["group"]].append(c)
    groups = {}
    for g in sorted(by_group):
        cs = by_group[g]
        if len(cs) < MIN_GROUP_QUALIFYING:
            continue
        etas = np.array([c["eta"] for c in cs])
        ci = boot_median_ci(etas, rng)
        prox = [c["eta_proximal"] for c in cs if c["eta_proximal"] is not None]
        dist = [c["eta_distal"] for c in cs if c["eta_distal"] is not None]
        groups[g] = {"archive": cs[0]["archive"], "cells": len(cs),
                     "median_eta": round(float(np.median(etas)), 4), "ci95": [round(ci[0], 3), round(ci[1], 3)],
                     "label": label(ci),
                     "cells_at_grid_bound": int(np.sum((etas <= ETA_GRID[0]) | (etas >= ETA_GRID[-1]))),
                     "secondary": {"median_eta_order_1_2": float(np.median(prox)) if prox else None,
                                   "median_eta_order_3plus": float(np.median(dist)) if dist else None,
                                   "median_fraction_within_10pct": float(np.median([c["within10"] for c in cs]))}}
    result = {"variant": variant, "cells_analysed": len(per_cell), "qualifying_groups": len(groups),
              "groups": groups}
    if not groups:
        result.update(status="Open", reason="no qualifying groups")
        return result
    archives = sorted({v["archive"] for v in groups.values()})
    arch_meds = {a: [v["median_eta"] for v in groups.values() if v["archive"] == a] for a in archives}
    pooled = float(np.median([v["median_eta"] for v in groups.values()]))
    boots = []
    for _ in range(N_BOOT):
        pick = rng.integers(0, len(archives), size=len(archives))
        boots.append(np.median(np.concatenate([arch_meds[archives[i]] for i in pick])))
    pooled_ci = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))
    counts = Counter(v["label"] for v in groups.values())
    n = len(groups)
    pooled_label = label(pooled_ci)
    if n < MIN_GROUPS_FOR_VERDICT:
        status = "Open"
    elif pooled_label == "consistent" and counts["consistent"] >= 2 * n / 3:
        status = "Supported"
    elif pooled_label == "deviating" and counts["deviating"] >= 2 * n / 3:
        status = "Refuted"
    elif counts["consistent"] >= 1:
        status = "Partial"
    else:
        status = "Open"
    result.update(pooled_median_eta=round(pooled, 3), pooled_ci95=[round(pooled_ci[0], 3), round(pooled_ci[1], 3)],
                  pooled_label=pooled_label, archives=len(archives), group_labels=dict(counts), status=status)
    return result


def full_run(offline: bool) -> tuple[dict, int]:
    neurons, info = load_metadata(offline)
    eligible = [n for n in neurons if eligible_by_metadata(n)]
    groups = defaultdict(list)
    for n in eligible:
        groups[group_key(n)].append(n)
    rng = np.random.default_rng(SEED)
    drawn = []
    for g in sorted(groups):
        members = sorted(groups[g], key=lambda n: n["neuron_id"])
        if len(members) < MIN_GROUP_ELIGIBLE:
            continue
        k = min(MAX_PER_GROUP, len(members))
        for i in sorted(rng.choice(len(members), size=k, replace=False)):
            drawn.append((g, members[int(i)]))
    exclusions = Counter()
    cells, manifest = [], []
    for i, (g, n) in enumerate(drawn):
        if i % 100 == 0:
            print(f"cell {i}/{len(drawn)}", file=sys.stderr)
        text, err = load_swc(n, offline)
        if text is None:
            exclusions[err] += 1
            continue
        triples, orders, why = cell_triples(text)
        if why is not None:
            exclusions[why] += 1
            continue
        cells.append({"group": g, "archive": n["archive"], "triples": triples, "orders": orders})
        manifest.append({"neuron_id": n["neuron_id"], "neuron_name": n["neuron_name"], "archive": n["archive"],
                         "reference_pmid": n.get("reference_pmid"),
                         "sha256": hashlib.sha256(text.encode()).hexdigest()})
    analyses = {v: analyse(cells, v, rng) for v in ANALYSES}
    primary = analyses["primary_5um"]["status"]
    flipped = [v for v in analyses if analyses[v]["status"] != primary]
    status = "Contested" if flipped else primary
    (HERE / "manifest.json").write_text(json.dumps(
        {"source": "NeuroMorpho.Org, CC BY 4.0", **info, "neurons": manifest}, indent=1))
    out = {"metadata": info, "eligible_by_metadata": len(eligible),
           "eligible_groups": len({g for g, _ in drawn}), "drawn": len(drawn),
           "analysed": len(cells), "exclusions": dict(exclusions),
           "analyses": analyses, "robustness_flips": flipped, "status": status}
    return out, 0 if status == "Supported" else 1


# ---- Amendment A1: secondary test on Blue Brain Open Data (does not set the status) ----

BB_BUCKET = "https://openbluebrain.s3.us-west-2.amazonaws.com/"
BB_PREFIX = "Experimental_Data/Reconstructed_morphologies/Categorized/Neurons/Rat/Cortex/SSCX/"
BB_DATA = DATA / "bluebrain"


def bb_list(offline: bool) -> list[tuple[str, int]]:
    """All (key, size) pairs under BB_PREFIX, from the S3 ListObjectsV2 API (cached)."""
    cache = BB_DATA / "listing.json"
    if cache.exists():
        return [tuple(x) for x in json.loads(cache.read_text())]
    if offline:
        raise SystemExit("no cached Blue Brain listing in data/bluebrain/; run without --offline first")
    import xml.etree.ElementTree as ET
    ns = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
    keys, token = [], None
    while True:
        q = {"list-type": "2", "prefix": BB_PREFIX}
        if token:
            q["continuation-token"] = token
        root = ET.fromstring(fetch(BB_BUCKET + "?" + urllib.parse.urlencode(q)))
        for c in root.findall("s3:Contents", ns):
            keys.append((c.find("s3:Key", ns).text, int(c.find("s3:Size", ns).text)))
        if root.findtext("s3:IsTruncated", namespaces=ns) != "true":
            break
        token = root.findtext("s3:NextContinuationToken", namespaces=ns)
    BB_DATA.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(sorted(keys)))
    return sorted(keys)


def bb_load(key: str, offline: bool) -> tuple[str | None, str | None]:
    path = BB_DATA / "swc" / hashlib.sha256(key.encode()).hexdigest()[:24]
    if path.exists():
        return path.read_text(errors="replace"), None
    if offline:
        return None, "not cached"
    try:
        raw = fetch(BB_BUCKET + urllib.parse.quote(key, safe="/"))
    except Exception as e:  # noqa: BLE001 - recorded as an exclusion reason
        return None, f"download failed: {type(e).__name__}"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)
    return raw.decode(errors="replace"), None


def bluebrain_run(offline: bool) -> dict:
    listing = bb_list(offline)
    by_group: dict[str, list[str]] = defaultdict(list)
    for key, _ in listing:
        rest = key[len(BB_PREFIX):].split("/")
        if len(rest) == 3 and rest[1] == "downloaded_from_nexus" and rest[2].lower().endswith(".swc"):
            by_group[rest[0]].append(key)
    name_count = Counter(k.rsplit("/", 1)[-1] for ks in by_group.values() for k in ks)
    dupes = {n for n, c in name_count.items() if c > 1}
    sizes = dict(listing)
    rng = np.random.default_rng(SEED)
    drawn = []
    for g in sorted(by_group):
        members = sorted(k for k in by_group[g] if k.rsplit("/", 1)[-1] not in dupes)
        if len(members) < MIN_GROUP_ELIGIBLE:
            continue
        k = min(MAX_PER_GROUP, len(members))
        for i in sorted(rng.choice(len(members), size=k, replace=False)):
            drawn.append((g, members[int(i)]))
    exclusions, cells, manifest = Counter(), [], []
    for i, (g, key) in enumerate(drawn):
        if i % 100 == 0:
            print(f"blue brain cell {i}/{len(drawn)}", file=sys.stderr)
        text, err = bb_load(key, offline)
        if text is None:
            exclusions[err] += 1
            continue
        triples, orders, why = cell_triples(text)
        if why is not None:
            exclusions[why] += 1
            continue
        # One institute: each m-type is its own resampling unit for the pooled interval.
        cells.append({"group": g, "archive": g, "triples": triples, "orders": orders})
        manifest.append({"key": key, "size": sizes.get(key), "m_type": g,
                         "sha256": hashlib.sha256(text.encode()).hexdigest()})
    analyses = {v: analyse(cells, v, rng) for v in ANALYSES}
    primary = analyses["primary_5um"]["status"]
    flipped = [v for v in analyses if analyses[v]["status"] != primary]
    (HERE / "manifest-bluebrain.json").write_text(json.dumps(
        {"source": "Blue Brain Open Data (s3://openbluebrain), CC BY 4.0", "prefix": BB_PREFIX,
         "files": manifest}, indent=1))
    return {"amendment": "A1", "role": "secondary; does not set the status",
            "source": BB_BUCKET + BB_PREFIX, "swc_files_listed": sum(len(v) for v in by_group.values()),
            "m_types_listed": len(by_group), "duplicate_names_dropped": len(dupes),
            "eligible_groups": len({g for g, _ in drawn}), "drawn": len(drawn), "analysed": len(cells),
            "exclusions": dict(exclusions), "analyses": analyses, "robustness_flips": flipped,
            "bluebrain_secondary_label": "Contested" if flipped else primary}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--bluebrain", action="store_true",
                        help="secondary test on Blue Brain Open Data only (Amendment A1)")
    args = parser.parse_args()
    t0 = time.time()
    validation = validate_estimator(n_cells=6 if args.quick else 20, seed=SEED)
    print("estimator validation:", json.dumps(validation))
    if args.quick:
        return 0 if validation["pass"] else 1
    if not validation["pass"]:
        print("estimator validation failed; not running the reanalysis")
        return 1
    if args.bluebrain:
        out = bluebrain_run(args.offline)
        record = {"result": "rall-3-2-rule", "mode": "bluebrain-secondary", "seed": SEED, "rall_eta": RALL_ETA,
                  "margin": list(MARGIN), "estimator_validation": validation, **out,
                  "python": platform.python_version(), "numpy": np.__version__,
                  "runtime_seconds": round(time.time() - t0, 1)}
        (HERE / "check-output-bluebrain.json").write_text(json.dumps(record, indent=2))
        p = out["analyses"]["primary_5um"]
        print(f"Blue Brain secondary label: {out['bluebrain_secondary_label']}; pooled eta "
              f"{p.get('pooled_median_eta')} CI {p.get('pooled_ci95')} over {p.get('qualifying_groups')} m-types")
        return 0
    out, code = full_run(args.offline)
    record = {"result": "rall-3-2-rule", "mode": "full", "seed": SEED, "rall_eta": RALL_ETA,
              "margin": list(MARGIN), "estimator_validation": validation, **out,
              "python": platform.python_version(), "numpy": np.__version__,
              "runtime_seconds": round(time.time() - t0, 1)}
    (HERE / "check-output.json").write_text(json.dumps(record, indent=2))
    p = out["analyses"]["primary_5um"]
    print(f"status: {out['status']}; pooled eta {p.get('pooled_median_eta')} "
          f"CI {p.get('pooled_ci95')} over {p.get('qualifying_groups')} groups")
    return code


if __name__ == "__main__":
    sys.exit(main())
