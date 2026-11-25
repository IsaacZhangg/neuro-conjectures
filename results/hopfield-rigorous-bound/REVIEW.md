# Review

| Date | Reviewer | Person or AI agent | Scope | Outcome |
| :-- | :-- | :-- | :-- | :-- |
| 2026-11-25 | Claude (drafted the proof and the check) | AI agent | Wrote PROOF.md; ran check.py (quick and full); certified A(τ₀) ∈ [0.1294899826213, 0.1294899826214]; arb cross-check agrees | Not an independent review |

Independent review: pending. A reviewer should check in particular:

- Lemma 4: the Sudakov–Fernique increment computation and the passage from
  finite index sets to the sphere.
- §5: the slack argument and that every error term is o(1) uniformly.
- That the cited theorems (Sudakov–Fernique, Gaussian concentration,
  Laurent–Massart, Talagrand's convex concentration) are stated in the form
  used; their theorem numbers were not re-checked against the books.
- check.py: that each series' ratio is non-increasing past the stopping
  point, which the tail bound relies on.
