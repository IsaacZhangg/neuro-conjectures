# Claim

## Statement tested

Target (open): raise the best rigorous lower bound on the storage capacity of
the Hebbian Hopfield model with independent unbiased ±1 patterns above 0.113,
toward the replica prediction α_c ≈ 0.138.

"Capacity" here is the retrieval notion used by the rigorous literature: with
probability tending to one, there is a stable fixed point of zero-temperature
single-spin-flip dynamics within Hamming distance δN of a stored pattern, for
some δ < 1/2. The strongest version (Newman's "energy barrier", called NLT by
Stojnic 2024) asks that every configuration at Hamming distance exactly
⌈δN⌉ from the pattern has strictly higher energy than the pattern; this
implies a stable fixed point inside the ball (Lemma 2 in PROOF.md).

What this result proves (our precise reading, stated in full in PROOF.md):

- **Theorem 1 (proved).** In the Hebbian model with one ±1 target pattern
  and M = αN − 1 *Gaussian* interference patterns, the energy barrier holds
  with probability at least 1 − e^{−cN} whenever
  α < A(τ) := ((δ(1−δ) − φ(τ)²)/φ(τ))² with δ = 1 − Φ(τ), for any τ > 0.
  At τ = 1.8259, A(τ) ≥ 0.12948 (certified with exact rational interval
  arithmetic). The supremum over τ is 0.1294900 (floating point), the same
  number as the first-level value in Stojnic (2024).
- **Theorem 2 (conditional).** The same bound holds for the standard model
  with ±1 patterns, for every pattern simultaneously, *if* one universality
  hypothesis (U) holds: the expected maximum of ‖Xσ‖ over the Hamming sphere
  is, to first order in N, the same for a ±1 matrix X as for a Gaussian one.
  (U) is not proved here or, as far as we found, anywhere.

## Where it comes from

- Feng, Shcherbina, Tirozzi (2001) state α ≥ 0.113 for ±1 patterns (abstract:
  "The probability of having a stable fixed point is one when α ≤ 0.113").
- Newman (1988) α ≥ 0.056, Loukianova (1994–97) α ≥ 0.071, Talagrand (1998)
  α ≥ 0.08 for the energy-barrier notion, as summarized in Stojnic (2024).
- Stojnic (2024) gives α_c(NLT,1) ≈ 0.129490 and α_c(AGS,1) ≈ 0.137906 for
  binary patterns, but works with Gaussian G and asserts the switch from
  binary to Gaussian "via concentrations and the Lindeberg principle" without
  a proof in the paper.

## Operational definitions

- Stable fixed point: σ with σ_i h_i(σ) > 0 for every i, where
  h_i(σ) = Σ_{j≠i} J_ij σ_j (for ±1 patterns, where energies can tie, ≥ 0
  with the convention that a spin with zero field keeps its value).
- Energy barrier at radius k: H(σ) > H(ξ) for every σ with d(σ, ξ) = k.
- Certified number: a rational lower bound obtained with exact rational
  arithmetic and explicit series remainders (check.py, no floating point in
  the certified path); cross-checked with Arb ball arithmetic when
  python-flint is installed.

## Amendment, 2026-11-25 (after the decision rule above was set)

Added Theorem 3 (PROOF-refstate.md). It keeps model G and the same Sudakov–Fernique bound, but compares the Hamming sphere to a lower-energy configuration σ₁ near ξ instead of ξ itself. This proves a stable fixed point near the target for α < 0.134716 (certified A₂(1.9422, 2.5312) ≥ 0.134716). This is a weaker notion than Newman's barrier against ξ, so it does not contradict Stojnic's NLT value 0.1295. It is still model G only. The ±1 status and the decision rule below are unchanged.

## Decision rule (set before writing the proof)

- **Proved** for Theorem 1 if PROOF.md gives a complete argument with every
  cited theorem stated, and check.py certifies A(τ₀) > 0.113.
- For the open ±1 target: **Partial** unless hypothesis (U) is proved. A
  bound for a different model, or a conditional bound, does not raise the
  rigorous ±1 bound, and the README must say so.
- **Disproved** is not available here: nothing in this folder can show that
  the ±1 capacity is below 0.113.

## Amendment, 2026-11-25 (novelty audit)

The claim and decision rule above are unchanged. A literature audit ([NOVELTY.md](NOVELTY.md)) found a preprint posted after this claim was written: Zhou, Hexagon 2611.00105 (2026-11-05), unrefereed. It claims ±1 fixed points up to α ≈ 0.1382. If it is correct, it answers this claim's question in the fixed-point sense by other means. What this folder adds that is still new is the energy barrier around every ±1 pattern at α < 0.12948 (Theorem 2, review pending). The previous rigorous barrier bound was 0.08.
