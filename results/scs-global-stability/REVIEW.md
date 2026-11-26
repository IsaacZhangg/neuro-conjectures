# Review

| Date | Reviewer | Person or AI agent | Scope | Outcome |
| :-- | :-- | :-- | :-- | :-- |
| 2026-11-25 | Claude (author) | AI agent | Proof written and self-checked; literature check by an AI research helper | Not an independent review |
| 2026-11-25 | Independent Claude referee (separate helper, did not write the proof) | AI agent | Every step of Theorems 1 and 2, constants, citations; numerical checks of the variance formula, edge mean, tail and Lipschitz constant | No fatal error. Fixed: boundary case and g ≥ 0 in Theorem 2(c); stated why the 1/√2 constant is sharp; diagonal-removal step stated as an inequality; global existence for general φ; notation for diagonal stability of J − I; citations (Bai–Yin, Geman, Tao, Lur'e–Postnikov). |
| 2026-11-26 | Independent Claude referee (separate helper, did not write BARRIER.md) | AI agent | BARRIER.md: Lemma 4, Theorem 3 Steps 1–4, corollary, remarks; numerical checks of the trace identity, the conditional variance of T_i and the Beta law of P_ii | Correct with minor fixes. Fixed: corollary strengthened to dV_D/dt > 0 near 0 by keeping half the margin; Step 1 constant renamed and described correctly; exponential rate in the no-self-coupling and elliptic remarks; scope of the Gibbs et al. open question stated (zero-mean Gaussian only); finite-N wording. |

No human review, at the repo owner's direction ("No human review.", 2026-11-25). The independent agent review above stands in for it. Under AGENTS.md an agent review cannot by itself earn ✅, so the mark stays 🟡 and the status stays Partial (proved for g < 1/√2, open for 1/√2 ≤ g < 1). The same holds for BARRIER.md.
