# Status note (2026-11-25)

**Update 22:05 UTC:** review resumed at the main chat's request; the result
is not to be posted or filed publicly until the review is done. Licenses now
checked (see SOURCES.md and DATA-LICENSE.md).

Earlier: parked. Isaac set a new direction for the collection: solve unsolved
conjectures rather than retest or verify existing ones. A new-scale test of
the Cuntz law counts as re-verification, so this result is not being taken
further for now.

## What was run

- Plan and decision rule frozen in commit 68a42ff, before any neuron was
  measured. Two recorded deviations in METHODS.md.
- Full run of `check.py` finished (442 s): FlyWire 783, 6,025 whole-neuron
  skeletons from 402 cell types; MICrONS v661, 1,309 dendritic skeletons from
  the V1 column. Output in check-output.json, per-neuron derived values in
  neuron-stats.csv, sources and SHA-256 in manifest.json.
- Result: β ≈ 0.53 in both datasets, below the 2/3 ± 0.1 window. FlyWire
  Refuted, MICrONS Contested (tips as n moves it inside). Status Contested.

## What was not done

- No independent review. An AI agent review was started and stopped
  unfinished when the target was parked; REVIEW.md is empty.
- No person has reviewed the code or numbers.
- No claim issue opened and no pull request; the branch is not merged.
- Not tried: the synapse form of the law, axons with presynaptic sites as n,
  hulls other than the convex hull, H01 human data, the hemibrain.

No raw data is committed. The derived per-neuron table holds only lengths,
counts and volumes computed by this check.
