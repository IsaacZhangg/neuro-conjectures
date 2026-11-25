# Collection design

This file sets the rules for what goes in the collection and how a result earns
its status. (It is about the collection, not about a website.)

## Scope

In scope: precise, falsifiable claims made in published neuroscience work, where
a reproducible check is possible with simulation, mathematics or public data.

- Theoretical and computational neuroscience: capacity bounds, dynamical
  regimes, optimal coding predictions, scaling laws derived from models.
- Quantitative empirical regularities: allometric and morphological scaling
  laws, population-code statistics, published effect sizes that public data
  can retest.

Out of scope:

- Claims that need new experiments, human subjects work, or data that is not
  openly licensed.
- Clinical claims, diagnosis or treatment advice.
- Vague claims with no stated quantity or test ("the brain is a predictive
  machine"). These can enter only after a specific, quantitative version is
  written down and attributed.

### Replication track

Some results retest claims that are already established (the first worked
example is one). They are labelled **replication** and counted separately from
new tests of open or contested claims. A replication must add something: an
independent implementation, a new dataset, a tighter scope, or a check of a
claim that has never been rerun.

## Claim kinds

| Kind | What counts as a check | Strongest status it can reach |
| :-- | :-- | :-- |
| Theorem | Complete written proof; Lean 4 formalization where practical | Proved / Disproved |
| Model claim | Deterministic simulation with fixed seeds, recorded output, stated finite-size limits | Supported / Refuted for the stated model |
| Empirical claim | Pre-registered reanalysis of named public datasets with exact versions | Supported / Refuted for the stated data |

A model claim says nothing about real brains unless a separate empirical result
links them. An empirical result says nothing beyond the species, preparation,
brain area and dataset it used.

## Status

Every result has exactly one status, written on its own line in its README as
`Status: <word>`.

- **Proved** / **Disproved**: theorems only. A complete argument has passed
  internal review.
- **Supported**: the check reproduces the claim within the stated scope and
  tolerance.
- **Refuted**: the check contradicts the claim within the stated scope, and the
  contradiction survives the robustness checks in METHODS.md.
- **Partial**: the claim was tested only in a narrower scope than stated (one
  dataset, finite sizes, one species, a special case).
- **Contested**: valid analyses disagree, or the result depends on an
  analysis choice that the original claim does not fix.
- **Open**: screened and claimed, no result yet. Open targets live in
  `docs/candidates.md`, not in the Results table.

Status words are never upgraded by rhetoric. "Supported in 7 of 7 sessions of
one dataset" is Partial unless the original claim was about that dataset.

## Screening a target

Before anyone works on a target, record in an issue:

1. The exact claim, quoted or paraphrased with its quantity, and the primary
   source with page, figure or equation number.
2. A literature check: has it already been retested, refuted or refined? Name
   what you searched.
3. Why it matters: what would change if it were supported or refuted.
4. The test: which model, proof strategy or dataset, and what result would
   count as support or refutation, decided in advance.
5. Data access: dataset, version or DOI, license.

## Pre-registration

For empirical claims, `METHODS.md` must be committed, with its decision rule,
before the analysis code reads the data. The git history is the record.
Deviations are allowed but must be listed under "Deviations from plan" with a
reason. This keeps reanalysis from turning into a search for a convenient
result.

## Review

A result moves from draft to the Results table only after:

- the check runs from a clean environment with pinned dependencies,
- someone other than the author has read the claim, sources and method and
  signed off in `REVIEW.md`,
- the README's scope statement matches what the check actually tested.

AI agents may draft, run and review. A review written by an agent says so, and
an agent review alone cannot earn ✅.

## Branches

`main` holds only reviewed results and accurate indexes. Work in progress lives
on `research` or short-lived branches. A publication pull request updates the
README table, CONTENTS.md and HISTORY.md together, and must pass
`python scripts/check_catalogue.py`.
