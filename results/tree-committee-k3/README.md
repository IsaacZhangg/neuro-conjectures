# Tree committee machine capacity, K = 3: a first proved lower bound

How many random patterns can a committee of three sign perceptrons, each
reading its own third of the input, store? The exact answer is open. The best
proved upper bound is 3.43 patterns per input weight, and before this result no
explicit proved lower bound was known. We prove that every load up to
**α = 1.55** is storable with probability tending to one, for Gaussian
patterns.

Status: Partial
Kind: Theorem
Reproducibility: 🟡 written proof and certified check with recorded output; not yet independently reviewed

## Result

**Theorem.** For the K = 3 tree committee machine with i.i.d. standard Gaussian
patterns and arbitrary ±1 labels, α_c(3) ≥ 31/20 = 1.55.

The proof is a three-stage construction. Branch 1 minimises a convex loss that
trades errors against margin. Branch 2 must be correct on every pattern that
branch 1 got wrong, and minimises the same kind of loss on the rest. Branch 3
only has to be correct on the patterns that branch 1 or branch 2 got wrong,
which is fewer than twice its dimension, so Cover's theorem finds it. The error
fractions of the first two branches come from the convex Gaussian min-max
theorem. Every numerical inequality is checked in Arb ball arithmetic; the
decisive one is that branch 3 must fit 1.99951 N < 2N patterns.

Gap left open: 1.55 ≤ α_c(3) ≤ 3.43. Replica theory (one step of replica
symmetry breaking) predicts about 3.0 to 3.14, and simulations suggest 2 to 2.75.

## Files

- [CLAIM.md](CLAIM.md): the open question and the exact statement proved.
- [PROOF.md](PROOF.md): the full proof.
- [SOURCES.md](SOURCES.md): references and searches.
- [METHODS.md](METHODS.md): how the construction and constants were chosen.
- [check.py](check.py): certified constants. `pip install -r requirements.txt`,
  then `python check.py` (well under a second).
- [check-output.json](check-output.json): recorded run.
- [sim_pipeline.py](sim_pipeline.py): finite-size simulation of the three
  stages (sanity check only; needs numpy, scipy and cvxpy).
- [REVIEW.md](REVIEW.md): review record.

## Narrower scope than the open question

- A lower bound, not the exact capacity. The open problem stays open.
- Gaussian patterns only. The usual ±1 inputs would need a universality step
  for the convex stages, which is not proved here.
- The construction is greedy and convex, so it cannot reach the true capacity.
  An idealised non-convex version (Gardner–Derrida margin laws, not rigorous)
  would give about 2.41.
- Agent-drafted proof (Claude), not yet independently reviewed.
