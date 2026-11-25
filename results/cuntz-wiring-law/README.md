# Cuntz wiring law on EM connectomes

Does total dendrite length scale as V^(1/3) · bp^(2/3), the wiring-optimal
scaling law of Cuntz, Mathy and Häusser (2012), in whole-cell EM
reconstructions from the fly and mouse connectomes?

Status: Open
Kind: Empirical claim
Reproducibility: ⬜ analysis plan only

## Result

Not run yet. The analysis plan and decision rule are frozen in
[CLAIM.md](CLAIM.md) and [METHODS.md](METHODS.md).

## Files

- [CLAIM.md](CLAIM.md): exact statement and decision rule.
- [SOURCES.md](SOURCES.md): primary source and later work.
- [METHODS.md](METHODS.md): data, analysis plan, robustness checks.
- [check.py](check.py): entry point. `python check.py` / `python check.py --quick`.
- check-output.json: recorded full run (added after the run).
- [REVIEW.md](REVIEW.md): review record.

## Scope and limits

One fly (FlyWire 783) and one mouse (MICrONS minnie65 v661). Branch-point form
of the law only; the synapse form was tested on EM data by Ferreira Castro and
Cardona. See METHODS.md for the known limits.
