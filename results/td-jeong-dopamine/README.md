# Can a TD model reproduce Jeong et al. 2022's dopamine results?

Open question: is there one temporal-difference reward prediction error model,
with a fixed state representation and one parameter set per animal, that passes
all eleven dopamine tests of Jeong et al. 2022 on the animals' own data, or can
one test be shown out of reach for a defined family of TD models?

Status: Open
Kind: Model claim
Reproducibility: ⬜ analysis plan only

## Result

Not run yet. CLAIM.md and METHODS.md were committed before any code read the
data.

## Files

- [CLAIM.md](CLAIM.md): exact statement and decision rule.
- [SOURCES.md](SOURCES.md): primary source and later work.
- [METHODS.md](METHODS.md): model family, tests, pass rules, search plan.
- [check.py](check.py): entry point. `python check.py` / `python check.py --quick`.
- check-output.json: recorded full run (added after the run).
- [REVIEW.md](REVIEW.md): review record.

## Scope and limits

Tests only the TD family defined in METHODS.md (linear TD(λ) over microstimulus
or serial-compound features of the event history, dopamine = TD error through a
fixed photometry kernel). A model that passes shows that the tests do not rule
out TD as such; it does not show that the brain runs that model. Licking is
taken from the data, not modelled.
