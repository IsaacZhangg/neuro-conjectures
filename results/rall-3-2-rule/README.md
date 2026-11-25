# Rall's 3/2 power rule at dendritic branch points

At a dendritic branch point, the parent diameter to the 3/2 power should equal
the sum of the daughter diameters to the 3/2 power; this tests whether real
dendrites in NeuroMorpho.org reconstructions do.

Status: Open
Kind: Empirical claim
Reproducibility: ⬜ analysis plan only

## Result

Not run yet. If the download fails with `CERTIFICATE_VERIFY_FAILED`, fix
the local certificate store (for example `pip install certifi` and
`export SSL_CERT_FILE=$(python -m certifi)`); do not turn verification off.

The full run has not been done. The plan and decision rule (CLAIM.md, METHODS.md) were committed
before any code read the data.

## Secondary test: Blue Brain rat cortex (does not set the status)

Planned in METHODS.md, Amendment A1, and committed before any Blue Brain file
was read. Data: 709 reconstructions drawn from 40 rat somatosensory cortex
cell types (m-types) in the Blue Brain Open Data bucket; 581 passed the file
filters, and 34 m-types had at least 10 passing cells.

- Pooled exponent: median 1.54, 95% interval 1.46 to 1.71 (bootstrap over
  m-types; one institute, so this interval is likely too narrow).
- Per cell type: 31 of 34 inconclusive, 3 deviating (L5 large and nest basket
  cells near 1.0 to 1.1, L6 bipolar cells near 2.8), none consistent. Cell
  exponents vary too much for 10 to 20 cells per type to settle a ±0.15
  margin.
- Label under the frozen rule: Open, and Contested across robustness checks:
  section-mean diameters give a pooled 1.78 (1.70 to 1.91) and 10 µm windows
  1.80 (1.73 to 1.89, Partial), while dropping identical triples leaves 1.53.
- The exponent rises away from the soma: in 30 of 34 m-types, bifurcations of
  order 3 and higher fit a larger exponent (median 1.72) than orders 1 and 2
  (median 1.40). This is a secondary analysis and was not tested for
  significance.
- Only about 16% of bifurcations per cell are within ±10% of Rall's rule.

Read this as: in this dataset the typical exponent sits near 3/2 on average,
but how near depends on how diameters are measured, and individual branch
points scatter widely. Recorded in check-output-bluebrain.json and
manifest-bluebrain.json; re-run with `python check.py --bluebrain` (about 6
minutes, downloads about 550 MB into `./data/`, ignored by git).
`--bluebrain --offline` reproduces the same output from the cache.

## Files

- [CLAIM.md](CLAIM.md): exact statement and decision rule.
- [SOURCES.md](SOURCES.md): primary source and later work.
- [METHODS.md](METHODS.md): data, inclusion rules, analysis plan, robustness checks.
- [check.py](check.py): entry point, added in the commit after the plan.
  `python check.py --quick` (under a second) validates the estimator on
  synthetic trees. `python check.py` downloads about 20 reconstructions per
  qualifying group from NeuroMorpho.org into `./data/` (ignored by git), so it
  takes hours and needs network access. `python check.py --offline` re-runs
  the analysis from that cache.
- manifest.json: neuron ids, archives, PubMed ids and SHA-256 of every analysed
  file (written by the full run).
- check-output.json: recorded full run (added after the run).
- check-output-bluebrain.json and manifest-bluebrain.json: recorded secondary
  run on Blue Brain Open Data (CC BY 4.0), with S3 key, size and SHA-256 of
  every analysed file.
- [REVIEW.md](REVIEW.md): review record.

## Scope and limits

Limited to NeuroMorpho.org reconstructions flagged as having diameters, from
control animals, in groups of at least 10 cells per archive, species and cell
type. Says nothing about cells imaged at resolutions where diameters are not
measured, or about axons.
