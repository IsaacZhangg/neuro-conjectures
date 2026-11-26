# Work in progress (paused 2026-11-26)

Not a result yet. Scratch scripts and the plan, saved so nothing is lost.

- bases.py: for k = 3, 4, 5 enumerates every full-support c = E^{-1}s
  (E an invertible k×k ±1 matrix, s a sign vector) and counts
  G_c = {ε ∈ {±1}^k : c·ε = ±1}. Found: k=3 all |G_c|/8 = 3/4; k=4 max 1/2;
  k=5 bases 5/8 (the ±1-coefficient c, 32 of them), 7/16, 5/16.
- inter.py r: pairwise intersections of the triple events and the 5-set
  ±1 events. Base 5/8 only for two triples sharing two indices with
  matching relative sign (48·C(r,4) such pairs); next bases 9/16.

Plan (CLAIM/METHODS not yet frozen):
1. Freeze CLAIM.md and METHODS.md for
   P(r,n) = 4C(r,3)(3/4)^n + (16C(r,5) − 48C(r,4))(5/8)^n + O(n^15 (9/16)^n),
   uniformly for 3 ≤ r ≤ n − 10n/ln n.
2. Proof: Bonferroni on triple and 5-set ±1 events; all other minimal
   supports by Odlyzko's pivot-row union bound with Erdős's two-antichain
   bound ρ_k ≤ (two largest binomials)/2^k (35/64 for k = 6, 7; ≤ 63/128 for
   k ≥ 8). Still to check: triple-intersection bases for the two
   configurations where all three pairs are 5/8-type (common pair, or three
   triples inside one 4-set) must be ≤ 9/16; rigorous tail sum for k ≥ 6.
3. check.py: exact r = 4 Möbius computation over all 2^16 row-pattern sets;
   expect coefficient 16 at (12/16)^n and −48 at (10/16)^n.
4. Independent agent review, catalogue checks, merge.
