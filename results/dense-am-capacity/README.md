# All-pattern capacity of dense associative memory with Gaussian patterns

Can a dense associative memory (modern Hopfield network) with e^{αN} Gaussian
patterns retrieve every pattern at once exactly up to the union-bound threshold
of Lucibello and Mézard, which is (log 2)/4 at large λ?

Status: Open
Kind: Theorem
Reproducibility: ⬜ analysis plan only

## Result

Not yet. Claim and proof plan frozen in CLAIM.md and METHODS.md.

## Files

- [CLAIM.md](CLAIM.md): exact statement and decision rule.
- [SOURCES.md](SOURCES.md): primary source and later work.
- [METHODS.md](METHODS.md): proof plan, constants, numerical support.
- [check.py](check.py): entry point. `python check.py` / `python check.py --quick`.
- check-output.json: recorded full run (added after the run).
- [REVIEW.md](REVIEW.md): review record.

## Scope and limits

Static signal-dominance definition of retrieval (Lucibello and Mézard eqs. 6 and
86), iid standard Gaussian patterns, N → ∞ with P = ⌊e^{αN}⌋. Not the dynamical
definition, not spherical or binary patterns.
