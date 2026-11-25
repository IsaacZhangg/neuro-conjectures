# Claim

## Statement tested

**Open question.** Let α_c(3) be the storage capacity of the tree committee
machine with K = 3 hidden sign units and non-overlapping receptive fields
(n = 3N inputs, output = majority of the three hidden signs), for m = αn random
patterns with random labels, as N → ∞. Find α_c(3). Any strictly tighter proved
bound on either side is progress.

**Statement proved here.** For i.i.d. N(0,1) pattern entries and labels that are
fixed or independent of the patterns, P(some weight vector stores all ⌊αn⌋
patterns) → 1 for every α ≤ 31/20. Hence α_c(3) ≥ 1.55 in that model.

## Where it comes from

- Upper bounds: Mitchison and Durbin (1989), about 5.42 for K = 3; Stojnic,
  arXiv:2312.08244 (random duality theory, 4.02, equal to the replica symmetric
  value); Stojnic, arXiv:2312.08257 (partially lifted random duality theory,
  3.43, the best proved upper bound we found).
- Non-rigorous: replica symmetric 4.02 (Barkai, Hansel and Sompolinsky 1992;
  Engel et al. 1992); one-step replica symmetry breaking about 3.0 to 3.14
  (Malzahn and Engel 1999 and references); simulations about 2 to 2.75.
- Lower side: we found no explicit proved lower bound at the scale of the
  phase transition (searches in SOURCES.md). A Cover-only construction gives
  4/3 (PROOF.md, remark in section 0); we found it nowhere in print either.

## Operational definitions

- "Stores": σ(ξ^μ) = y_μ for every μ, with σ = sign of the sum of the three
  hidden signs.
- "With probability tending to one": over the Gaussian patterns, for every
  fixed label sequence.
- Capacity lower bound: α_c(3) ≥ α means storability at load α with
  probability → 1.

## Decision rule

A theorem has no data to peek at, so the rule is about the argument.

- Proved (as a bound) if PROOF.md gives a complete argument and every
  numerical inequality it uses is certified by check.py (interval arithmetic,
  exit code 0).
- The result folder is labelled **Partial** whatever happens, because it bounds
  the capacity instead of determining it, and covers Gaussian patterns only.
- Refuted (withdrawn) if review finds a gap that cannot be closed.
