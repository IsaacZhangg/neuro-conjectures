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
- [REVIEW.md](REVIEW.md): review record.

## Scope and limits

Limited to NeuroMorpho.org reconstructions flagged as having diameters, from
control animals, in groups of at least 10 cells per archive, species and cell
type. Says nothing about cells imaged at resolutions where diameters are not
measured, or about axons.
