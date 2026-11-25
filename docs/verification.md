# How results are checked

Every result folder must let a stranger rerun the check and see the same
numbers. This page lists what each kind of check must contain.

## Common to every result

- `README.md` with the claim in one sentence, `Status:` line, scope, and limits.
- `CLAIM.md`: the exact statement being tested and the decision rule.
- `SOURCES.md`: primary source with page, equation or figure number, plus any
  later papers that retest, refine or dispute it, and the date sources were
  checked.
- `METHODS.md`: model or dataset, parameters, analysis plan, robustness checks,
  deviations from plan.
- `check.py` (or another single entry point) that exits 0 when the decision
  rule passes and 1 when it fails, and writes `check-output.json`.
- `check-output.json`: the recorded output of the full run, with software
  versions, seed and runtime.
- `REVIEW.md`: who reviewed what, and whether the reviewer was a person or an
  AI agent.
- Pinned dependencies (`requirements.txt` at the root, or one in the folder).

A `--quick` mode that finishes in under a minute is required so CI can run
every result on every pull request. Quick mode is a smoke test that the code
runs; only the full run sets the status.

## Theorems

- A complete written proof in `PROOF.md`, stating every assumption.
- Optional `lean/` formalization. Mark clearly what the Lean code proves, what
  it assumes as an axiom, and where it differs from the paper's statement.
- Numerical spot checks are welcome but never replace the proof.

## Model claims (simulation)

- Fixed random seeds, deterministic code, recorded output.
- State the exact model: equations, update rule, parameters, initial
  conditions. Small modelling choices (synchronous vs asynchronous updates,
  self-couplings) can change results, so name them.
- Finite-size and finite-time limits. If the claim is asymptotic (N → ∞), show
  several sizes and say how you extrapolated and why.
- The decision rule (target value and tolerance) is set in `CLAIM.md` before the
  full run.

## Empirical claims (reanalysis of public data)

- Name each dataset with a persistent identifier (DOI, DANDI ID, OpenNeuro
  accession, NeuroMorpho archive) and the exact version or download date.
- `METHODS.md` and the decision rule are committed before the analysis reads the
  data. The git history is the evidence.
- Separate the planned analysis from exploratory follow-ups. Only the planned
  one sets the status.
- Report effect sizes with uncertainty (confidence or credible intervals,
  bootstrap over animals or sessions, not over neurons when neurons share an
  animal).
- Robustness: at least one reasonable alternative analysis choice (a different
  preprocessing step, fitting range or inclusion rule). If the conclusion flips,
  the status is Contested.
- Do not commit raw data. Commit a download script with checksums.
- Respect dataset licenses and terms of use. Never use data that identifies a
  person.

## What a check cannot show

- A simulation shows what a model does, not what a brain does.
- A reanalysis shows what one dataset shows. It cannot establish mechanism or
  causation on its own.
- A failed replication may reflect a difference in data or methods rather than
  a false claim. Say which differences you know about.
