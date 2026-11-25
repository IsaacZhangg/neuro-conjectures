# Rigorous Hopfield capacity bound, Gaussian-interference model

Can it be proved that a Hebbian Hopfield network keeps a stable state near
a stored pattern for loads above 0.113, the best rigorous bound on record?

Status: Partial
Kind: Theorem
Reproducibility: 🟡 written proof and certified check with recorded output; not yet independently reviewed

## Result

**Proved (Theorem 1), for a narrower model.** With one ±1 target pattern and
Gaussian interference patterns, every configuration at Hamming distance
⌈0.0339 N⌉ from the target has higher energy than the target, and a stable
fixed point with overlap at least 0.932 exists, with probability
1 − e^{−cN}, for every load **α < 0.12948**. The constant is certified with
exact rational arithmetic: A(1.8259) ∈ [0.1294899826213, 0.1294899826214].
The proof uses the Sudakov–Fernique comparison plus Gaussian concentration
([PROOF.md](PROOF.md) §1–5). The value equals Stojnic's (2024) first-level
NLT capacity 0.129490; what is new here is a complete written proof of it
for this model, not a new constant.

**Not proved: the standard ±1 model.** This does **not** raise the
rigorous bound for the usual Hopfield model with ±1 patterns, which stays at
0.113 (Feng, Shcherbina, Tirozzi 2001) as far as we found. Theorem 2 gives
the same 0.12948 for ±1 patterns, for every pattern at once, *only if* a
universality hypothesis (U) holds: the expected maximum of ‖Xσ‖ over the
Hamming sphere is the same to first order for ±1 and Gaussian X. PROOF.md §6
explains why the standard ±1-to-Gaussian tools do not give (U).

## Files

- [CLAIM.md](CLAIM.md): the statements and decision rule.
- [PROOF.md](PROOF.md): complete proof of Theorem 1; conditional Theorem 2;
  the gap.
- [check.py](check.py): certified computation of the constant (standard
  library only) plus optional Arb cross-check and exploratory simulation.
  `python check.py` or `python check.py --quick` (both under a second
  without the simulation; a few seconds with numpy).
- [check-output.json](check-output.json): recorded full run.
- [SOURCES.md](SOURCES.md), [METHODS.md](METHODS.md),
  [REVIEW.md](REVIEW.md) (independent review pending).

## Scope and limits

- Theorem 1 is about Gaussian interference patterns, not ±1 ones. It is a
  special case relative to the target, hence Partial.
- The proof is agent-drafted and has not been independently reviewed. Per
  DESIGN.md, "Proved" needs internal review, so even Theorem 1 is not
  labelled Proved yet.
- The energy-barrier notion is stronger than the AGS basin notion behind
  0.138; the first-level NLT value 0.1295 is not expected to reach 0.138.
- The theorem numbers of cited textbook results in PROOF.md §7 were not
  re-checked against the books.
- About the model, not about the brain.
