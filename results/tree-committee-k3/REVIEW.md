# Review

| Date | Reviewer | Person or AI agent | Scope | Outcome |
| :-- | :-- | :-- | :-- | :-- |
| 2026-11-25 | Claude (author) | AI agent | Self-check of PROOF.md against check.py output; simulation sanity check | No gap found by the author. Not an independent review. |
| 2026-11-26 | Claude (author) | AI agent | Re-computed E[D²] and E[(G+D)D] for both stages at the ends of the certified root brackets by direct numerical quadrature of the piecewise displacements (scipy), independent of the closed forms in check.py | Agrees with the Arb values to at least 6 digits: stage 1 P = 1.905e-4 to 1.911e-4, stage 2 P = 1.2769e-3 to 1.2773e-3, stage-3 load 1.999505. Not an independent review. |

Independent review pending. Points a reviewer should check first:

1. Lemma 2, Step 2: the exact-penalty step and the bound on the multipliers.
2. Lemma 2, Step 4: uniform convergence of the auxiliary problem and the sign of
   its derivative in r at r = 1.
3. Lemma 2, Step 5: the deviation argument near r = 1.
4. Independence: H, R depend only on branch-1 data; T only on branch-1 and
   branch-2 data.
