# All-pattern capacity of dense associative memory with Gaussian patterns

Can a dense associative memory (modern Hopfield network) with e^{αN} Gaussian
patterns retrieve every pattern at once exactly up to the union-bound threshold
of Lucibello and Mézard, which is (log 2)/4 at large λ?

Status: Proved
Kind: Theorem
Reproducibility: 🟡 runnable check with recorded output, agent-reviewed, not yet run by CI

## Result

Yes, under the static retrieval criterion. For every λ ∈ (0, ∞], the
all-pattern threshold equals the union bound: α_c(λ) = min(α_a(λ), (log 2)/4),
with a sharp transition (P(all retrieved) → 1 below, → 0 above). So
α_c = (log 2)/4 for all λ ≥ λ_x ≈ 0.70092, including λ = ∞, and
α_c(λ) = α_a(λ) < (log 2)/4 below λ_x (α_a is given by a Lambert-W formula).
The matching upper bound comes from one competitor beating a short pattern of
squared norm about N/2, the same configuration the union bound picks, found
by a second-moment estimate. Proof in [PROOF.md](PROOF.md).

## Files

- [CLAIM.md](CLAIM.md): exact statement and decision rule.
- [SOURCES.md](SOURCES.md): primary source and later work.
- [NOVELTY.md](NOVELTY.md): what was already known, what is new, verdict.
- [METHODS.md](METHODS.md): proof plan, constants, numerical support.
- [check.py](check.py): entry point. `python check.py` / `python check.py --quick`.
- [PROOF.md](PROOF.md): the proof.
- [check-output.json](check-output.json): recorded full run.
- [REVIEW.md](REVIEW.md): review record.

## Scope and limits

Static signal-dominance definition of retrieval (Lucibello and Mézard eqs. 6 and
86), iid standard Gaussian patterns, N → ∞ with P = ⌊e^{αN}⌋. Not the dynamical
definition (convergence of the update rule from a nearby start, which is how
the paper first defines retrieval), not spherical or binary patterns. No Lean
formalization; two independent reviews, both by AI agents, and no human
review.
