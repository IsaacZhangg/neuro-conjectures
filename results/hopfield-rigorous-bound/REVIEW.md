# Review

| Date | Reviewer | Person or AI agent | Scope | Outcome |
| :-- | :-- | :-- | :-- | :-- |
| 2026-11-25 | Claude (drafted the proof and the check) | AI agent | Wrote PROOF.md; ran check.py (quick and full); certified A(τ₀) ∈ [0.1294899826213, 0.1294899826214]; arb cross-check agrees | Not an independent review |
| 2026-11-25 | Claude (drafted Theorem 3 in a helper session; checked by the main session) | AI agent | Read PROOF-refstate.md step by step (decomposition, exact identity, tail bounds, slack constants); reran check_refstate.py with arb | No gap found; not an independent review |
| 2026-11-26 | Claude (separate reviewer agent, did not write the proof) | AI agent | Independent adversarial review of PROOF-pm1.md §2–§7 against PROOF.md; reran check_pm1.py (20/20 pass, arb agrees); recomputed δ₀, φ₀, A(τ₀), sup A, A_H(δ₀) with its own script; grid check of Proposition 4.4 | Passes. No fatal or fixable gap. Six cosmetic issues, fixed on 2026-11-26 (see below) |

### Independent agent review of PROOF-pm1.md (2026-11-26)

Verdict: passes. Points checked and found sound:

- Step 4: the increment identity and the sign of the Sudakov–Fernique
  interpolation bracket; the common mean is deterministic once X_T is fixed.
- Step 3: Lindeberg is applied only to the smooth soft maximum; its
  third-derivative bound uses only the caps on the swapped rows; the error is
  o(N) with N → ∞ first.
- The net enters only through the soft maximum, never a union bound. The two
  union bounds cost N·O(√(ε log(1/ε))), which vanishes as ε → 0, and the limits
  are taken in the right order (Corollary 6.1).
- The union over patterns needs no independence.
- Cases A and B together cover all of (0, 0.12948).
- The zero-diagonal convention is consistent with the barrier criterion (★).
- check_pm1.py certifies exactly the inequalities the proof uses, with outward
  rounding.

Cosmetic issues raised and fixed: the list of quoted results in §0 (the net
bound is proved in Lemma 3.3; Taylor's theorem and the Mills bound were
missing); reused symbols (ψ, H, Y, Φ), now listed in §0; an overstated remark in
§7 on the Hoeffding route; PROOF.md §6 still described (U) as open, now with a
pointer to PROOF-pm1.md. Still not re-checked: the theorem numbers of the
quoted textbook results, and the arXiv numbers cited in PROOF.md §6.

No person has reviewed the proofs yet.

### Earlier checklist for PROOF.md (Theorem 1)

A reviewer should check in particular:

- Lemma 4: the Sudakov–Fernique increment computation and the passage from
  finite index sets to the sphere.
- §5: the slack argument and that every error term is o(1) uniformly.
- That the cited theorems (Sudakov–Fernique, Gaussian concentration,
  Laurent–Massart, Talagrand's convex concentration) are stated in the form
  used; their theorem numbers were not re-checked against the books.
- check.py: that each series' ratio is non-increasing past the stopping
  point, which the tail bound relies on.
