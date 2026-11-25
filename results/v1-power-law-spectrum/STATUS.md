# Status: parked (2026-11-25)

Parked before any real data were analysed. The project's direction changed to
working on unsolved conjectures rather than re-testing published findings, and
this target is a re-test of Stringer et al. 2019.

## What was done

- Screened the claim and confirmed sources (SOURCES.md).
- Froze CLAIM.md and METHODS.md (decision rule, bands, forward model) before
  any code read neural data.
- Wrote check.py; `python check.py --quick` (synthetic only) passes.
- Downloaded the planned Allen Brain Observatory sessions (456 two-photon
  session B experiments, 32 Neuropixels sessions) to a cloud cache only.
  Nothing from them is committed, and no spectrum was computed from them.
- Synthetic validation of the frozen method found problems (broken power law
  wrongly preferred, bootstrap intervals missing the true α), which led to the
  pre-data amendment recorded under "Deviations from plan" in METHODS.md.
- The synthetic generator was replaced by one whose population spectrum is
  exactly known. With the amendment, partial re-validation (6 of 24 planned
  runs) gave:

| true α | K | neurons | α̂ | 95% interval | broken PL clearly better |
|---|---|---|---|---|---|
| 1.04 | 3000 | 4000 | 1.04 | [1.04, 1.08] | no |
| 1.04 | 3000 | 4000 | 1.04 | [0.94, 1.08] | no |
| 1.04 | 3000 | 4000 | 1.04 | [1.06, 1.08] | no |
| 1.04 | 10000 | 12000 | 1.04 | [1.04, 1.06] | no |
| 1.04 | 10000 | 12000 | 1.04 | [1.04, 1.04] | no |
| 1.04 | 10000 | 12000 | 1.04 | [1.04, 1.06] | no |

Point estimates are right and the false positives are gone, but the bootstrap
intervals are too narrow and sometimes miss the truth (third row), so the
interval method still needs work before any real-data run.

## Not done

- Remaining validation cases (α 1.4, 0.8, 1.2, a broken-power-law truth for
  power).
- A fix for the interval coverage.
- The full run on real data, check-output.json, REVIEW.md entries, catalogue
  status, and the public claim issue (a draft exists outside the repo; not
  opened).
