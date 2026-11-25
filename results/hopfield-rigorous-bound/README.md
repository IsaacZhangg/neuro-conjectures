# Rigorous Hopfield capacity bound, Gaussian-interference model

Can it be proved that a Hebbian Hopfield network keeps a stable state near
a stored pattern for loads above 0.113, the best rigorous bound on record?

Status: Partial
Kind: Theorem
Reproducibility: 🟡 written proofs and certified checks with recorded output; ±1 proof under independent review

## Result

**Novelty audit (2026-11-25), read first: [NOVELTY.md](NOVELTY.md).** A
preprint posted on 2026-11-05 by Zhou claims the zero-temperature AGS
capacity for ±1 patterns lies in [0.138198, 0.138201], with fixed points
near every pattern. It is unrefereed and AI-written. If it is correct, it
supersedes this folder's original goal of a ±1 fixed-point bound above
0.113, and it supersedes Theorem 3. The part that still looks new is the
firm energy barrier around each ±1 pattern (the "NLT" notion). The best
rigorous bound in print for that notion is 0.08 (Talagrand 1998).

**Theorem 2, ±1 patterns, energy barrier around every pattern: proof
written, review pending.** With ±1 patterns, with probability 1 − e^{−cN},
every pattern is surrounded by an energy barrier at Hamming radius
⌈0.0339 N⌉, so it has a stable fixed point with overlap at least 0.932.
This holds for every load **α < 0.12948**. [PROOF-pm1.md](PROOF-pm1.md)
proves the universality step (U) that PROOF.md §6 left open. It uses a
dissection of the maximizing pattern directions, a Lindeberg swap on the
capped rows only, a conditional Sudakov–Fernique comparison, and a scalar
Hoeffding and Lindeberg bound. Its constants are certified exactly by
[check_pm1.py](check_pm1.py). The constant equals Stojnic's (2024)
non-rigorous first-level NLT value 0.1294899. The previous rigorous NLT
bound was 0.08, and random duality theory predicts the true value is about
0.1298. No independent review has finished yet, so this is not labelled
Proved.

**Theorem 1, Gaussian interference: proved.** This is the same barrier
for one ±1 target with Gaussian interference patterns, for α < 0.12948.
The constant is certified: A(1.8259) ∈ [0.1294899826213, 0.1294899826214]
([PROOF.md](PROOF.md) §1–5).

**Theorem 3, Gaussian interference, fixed point: proved, superseded if
Zhou holds.** The proof compares the sphere with a lower-energy reference
state near the target. It gives a stable fixed point for every α <
0.134716, certified at A₂(1.9422, 2.5312) ∈ [0.134716551489,
0.134716551490] ([PROOF-refstate.md](PROOF-refstate.md)). This uses the
weaker fixed-point notion, which Zhou's claimed 0.1382 covers.

## Files

- [CLAIM.md](CLAIM.md): the statements and decision rule.
- [PROOF.md](PROOF.md): complete proof of Theorem 1; conditional Theorem 2;
  the gap.
- [PROOF-refstate.md](PROOF-refstate.md): Theorem 3, the reference-state refinement.
- [PROOF-pm1.md](PROOF-pm1.md): proof of (U), making Theorem 2 unconditional for ±1 patterns (review pending); [check_pm1.py](check_pm1.py) certifies its constants.
- [NOVELTY.md](NOVELTY.md): literature audit, what is known and what is open.
- [drafts/](drafts/): the unreviewed planning draft behind PROOF-pm1.md.
- [check.py](check.py) (with [check_refstate.py](check_refstate.py)): certified computation of both constants (standard
  library only) plus optional Arb cross-check and exploratory simulation.
  `python check.py` or `python check.py --quick` (both under a second
  without the simulation; a few seconds with numpy).
- [mc_refstate.py](mc_refstate.py): Monte Carlo sanity check of the Theorem 3 energy identity (exploratory).
- [check-output.json](check-output.json): recorded full run.
- [SOURCES.md](SOURCES.md), [METHODS.md](METHODS.md),
  [REVIEW.md](REVIEW.md) (independent review pending).

## Scope and limits

- Theorems 1 and 3 are about Gaussian interference patterns. Theorem 2 covers
  ±1 patterns but is still under review, and its constant 0.12948 sits below
  the predicted NLT capacity of about 0.1298, hence Partial.
- The proof is agent-drafted and has not been independently reviewed. Per
  DESIGN.md, "Proved" needs internal review, so even Theorem 1 is not
  labelled Proved yet.
- The energy-barrier notion is stronger than the AGS basin notion behind
  0.138; the first-level NLT value 0.1295 is not expected to reach 0.138.
- The theorem numbers of cited textbook results in PROOF.md §7 were not
  re-checked against the books.
- About the model, not about the brain.
