# Rall's 3/2 power rule at dendritic branch points

At a dendritic branch point, the parent diameter to the 3/2 power should equal
the sum of the daughter diameters to the 3/2 power; this tests whether real
dendrites in NeuroMorpho.org reconstructions do.

Status: Open
Kind: Empirical claim
Reproducibility: ⬜ analysis plan only

## Result

Not run yet. The plan and decision rule (CLAIM.md, METHODS.md) were committed
before any code read the data.

## Files

- [CLAIM.md](CLAIM.md): exact statement and decision rule.
- [SOURCES.md](SOURCES.md): primary source and later work.
- [METHODS.md](METHODS.md): data, inclusion rules, analysis plan, robustness checks.
- [check.py](check.py): entry point (added after the plan commit).
- check-output.json: recorded full run (added after the run).
- [REVIEW.md](REVIEW.md): review record.

## Scope and limits

Limited to NeuroMorpho.org reconstructions flagged as having diameters, from
control animals, in groups of at least 10 cells per archive, species and cell
type. Says nothing about cells imaged at resolutions where diameters are not
measured, or about axons.
