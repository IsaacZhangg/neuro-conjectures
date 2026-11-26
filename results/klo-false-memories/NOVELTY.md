# Novelty check

Checked 2026-11-26. Searches: Exa web search (the Exa plugin server failed to
connect, the Exa search tool worked), Firecrawl paper search, arXiv full text
of the papers below, and their citing papers from 2023 to now.

## The conjecture (exact statement and source)

Let v_1, …, v_r be independent and uniform in {−1,1}^n and let P(r,n) be the
probability that their real span contains a ±1 vector other than ±v_1, …,
±v_r. For the projection (pseudo-inverse) Hopfield rule, the stored patterns
are the ground states of H(x) = −½ xᵀP_V x on {−1,1}^n, and the other ground
states are exactly the ±1 vectors in the span ("false memories").

**Kalai–Linial–Odlyzko (KLO) conjecture**, stated in Odlyzko, J. Combin.
Theory A 47 (1988) 124–133, and extended to the endpoint by Kahn, Komlós and
Szemerédi (KKS), J. AMS 8 (1995) 223–240, Section 4: for all r ≤ n − 1,

    P(r,n) = (1 + o(1)) · 4·C(r,3)·(3/4)^n.

## Verdict

**Partly known.** The sharp asymptotic is proved for r ≤ n − C. The weak form
(P → 0) is proved for every r ≤ n − 1. Only the sharp asymptotic in the last
window n − C < r ≤ n − 1 is open. We judge that window out of reach here (see
below) and retarget at an open refinement inside the proved range: the exact
second-order term.

## What was already known

- **Odlyzko 1988**, Theorem: P(r,n) = 4·C(r,3)·(3/4)^n + O((7/10)^n),
  uniformly for r ≤ n − 10n/ln n. Prop. 2.3: the three-vector part equals
  4·C(r,3)·(3/4)^n + O(r^4 (5/8)^n). He writes that the (7/10)^n error "can
  be substantially improved with additional effort".
- **KKS 1995**, Corollary 4(b): the same leading term for r ≤ n − C, with
  error O((7/10)^n). They conjecture r ≤ n − 1 suffices.
- **Irmatov**, arXiv:2405.05082 (2024): claims r ≤ n − 1 with error
  O((5/8 + o(1))^n). The proof rests on his own claimed singularity
  asymptotics (arXiv:2004.03400). Xie (below) writes that "several issues in
  the argument remain unresolved". We treat it as unproved.
- **Xie**, arXiv:2609.29756 v2 (Nov 20 2026), Theorem 3: for all
  1 ≤ r ≤ n − 1, the span contains no other ±1 vector with probability
  1 − O(n^{−99}). This settles the threshold r = n − 1 (the weak KLO form)
  but not the sharp asymptotic. Theorem 7 proves the KKS dependence
  conjecture 2·C(r,2)·2^{−n} + O(2^{−n}e^{−cn}) for r ≤ n − 1.
- **Hu, Song and Wu**, arXiv:2609.39463 (Oct 2026): claim the Rademacher
  singularity asymptotic. Not about spans; not yet refereed.
- Related, different questions: Ferber, Jain and Zhao, CPC 2021 (Hadamard
  matrices, quotes KKS Cor. 4(b)); Hunter, Kwan, Sauermann and Sawhney, IMRN
  2026; Bruck and Roychowdhury 1990 (spurious memories of the Hebbian rule,
  a different model); Voigt and Ziegler, math/0308050 (0/1 vectors).

## What is open

1. The sharp KLO asymptotic for n − C < r ≤ n − 1. Xie's Remark 8 notes that
   this precision at r = n would give the singularity asymptotic; the window
   just below n is of the same kind of difficulty, and Xie's cokernel method
   gives only polynomial error. Not attempted here.
2. **The error term inside the proved range.** Nobody states the next term.
   Odlyzko proves O((7/10)^n); KKS O((7/10)^n); Irmatov's (5/8)^n order is
   only claimed. We found no paper that gives the coefficient of (5/8)^n.

## Target here (PUSH)

Prove, uniformly for 3 ≤ r ≤ n − 10n/ln n,

    P(r,n) = 4·C(r,3)·(3/4)^n + (16·C(r,5) − 48·C(r,4))·(5/8)^n
             + O(n^15 (9/16)^n).

This replaces Odlyzko's O((7/10)^n) error by an exact second-order term,
confirms that the (5/8)^n order in Irmatov's claim is the true order (it is
attained, with a coefficient that is negative for r < 19 and positive for
r > 19), and gives the false-memory rate of the projection rule to second
order. Solving counts as a complete proof plus exact checks (CLAIM.md).

## Not covered

- The window n − C < r ≤ n − 1 (item 1 above).
- The range n − 10n/ln n < r ≤ n − C proved by KKS; our proof uses
  Odlyzko's counting and does not reach it.
