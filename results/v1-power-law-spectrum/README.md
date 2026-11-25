# Power-law eigenspectrum of visual cortex responses (α ≈ 1)

The variance of the n-th principal component of mouse visual cortex responses
to natural images falls off as 1/n^α with α just above 1.

Status: Open
Kind: Empirical claim
Reproducibility: ⬜ analysis plan only

## Result

Not run yet. The analysis plan and decision rule were committed before any
code read the data.

## Files

- [CLAIM.md](CLAIM.md): exact statement, our reading, and decision rule.
- [SOURCES.md](SOURCES.md): primary source, later work, datasets, searches.
- [METHODS.md](METHODS.md): data, analysis plan, robustness checks.
- [check.py](check.py): entry point. `python check.py` / `python check.py --quick`.
- check-output.json: recorded full run (added after the run).
- [REVIEW.md](REVIEW.md): review record.

## Scope and limits

- Mouse visual cortex only (six areas), Allen Brain Observatory two-photon and
  Neuropixels recordings, 118 natural images.
- With 118 images the test can only ask whether the observed spectrum matches
  what a power-law population spectrum would give for 118 images, under the
  forward-model assumptions in METHODS.md. It cannot see ranks above 117,
  where most of the original fit range (11 to 500 of 2,800) lies.
- The original Stringer et al. data could not be downloaded from the
  environment this was run in, so this is not a reanalysis of those data.
