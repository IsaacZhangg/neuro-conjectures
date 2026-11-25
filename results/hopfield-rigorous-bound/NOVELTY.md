# Novelty audit

Literature search done on 2026-11-25. Tools: Exa web search, the Firecrawl research index (arXiv, PubMed, bioRxiv) from 2020 to now, and papers citing Stojnic (2024). arXiv's own site turned away direct requests, so arXiv papers were read through these indexes. Google Scholar was not reachable from here.

**Verdict: partly known.** Most of what this folder set out to do was done by a preprint posted on 2026-11-05. That preprint has not been refereed. One part of this folder still looks new: the rigorous bound for the firm energy barrier around each ±1 pattern.

## What this folder states

All statements use Hebbian couplings with M = αN patterns.

- **(T1)** Gaussian interference patterns: there is an energy barrier at Hamming radius ⌈0.0339N⌉ around the target, for every α < 0.12948. This is proved.
- **(T2)** ±1 patterns: the same barrier holds around every pattern at once, for every α < 0.12948. The proof is in [PROOF-pm1.md](PROOF-pm1.md). It has been written out and certified, but no independent review has finished yet.
- **(T3)** Gaussian interference patterns: there is a stable fixed point near the target, for every α < 0.134716. This is proved.

The barrier in T1 and T2 is the "NLT" notion (Newman, Loukianova, Talagrand). Every configuration on the sphere has strictly higher energy than the pattern itself. A barrier of this kind implies a fixed point inside the sphere. T3 and the 0.138 capacity use the weaker "AGS" notion, which asks only for an energy well near the pattern.

## What is known

### Older rigorous results

| Notion | Patterns | Bound | Source |
| :-- | :-- | :-- | :-- |
| NLT barrier | ±1 | α ≥ 0.056 | Newman (1988), Neural Networks 1:223 |
| NLT barrier | ±1 | α ≥ 0.071 | Loukianova (1997), PTRF 107:161 |
| NLT barrier | ±1 | α ≥ 0.08 | Talagrand (1998), PTRF 110:177 |
| Barrier around some point inside the ball | ±1 | α ≥ 0.113 | Feng, Shcherbina, Tirozzi (2001), CMP 216:139 |

Up to October 2026, every survey we found still names 0.08 as the best NLT bound. Examples are Stojnic, arXiv:2610.38437 (2026), and the Upgang thesis (2020).

### Non-rigorous predictions

Stojnic, arXiv:2403.01907 (2024), uses random duality theory. His first lifting level gives an NLT capacity of 0.1294899, the same constant as T1 and T2. His second level gives 0.12979. He also recovers the replica AGS values: 0.137906 (replica symmetric) and 0.138186 (one step of replica symmetry breaking).

### New since 2026, all unrefereed

- **Zhou, "The critical capacity of the Hopfield model".** Hexagon preprint 2611.00105v1, posted 2026-11-05. Code and certificates are at github.com/OctopusTide/hopfield-capacity. The proofs were written mainly by AI; the author checked parts. What it claims:
  - The zero-temperature AGS capacity for ±1 patterns lies in [0.138198129, 0.13820106].
  - Binary patterns and Gaussian rows have the same wells, at every temperature.
  - For the minimum energy at fixed overlap, a well holds around all patterns at once with probability 1 − e^{−cN} (its Theorem 1.2(e)).
  - It relies on Ali, Boursier, Chen, Ho, Issa and Mourrat, "Free energy of the Hopfield model" (arXiv:2611.03306, 2026-11-02). That paper identifies the free energy for Bernoulli patterns by reducing them to Gaussian ones.
- **Camilli, Contucci and Mingione**, arXiv:2610.32618: the Hopfield model with Gaussian patterns.

## What this means for each statement

- **T3 (fixed point for α < 0.134716, Gaussian interference): superseded, if Zhou is correct.**
  - A well in the minimum energy at fixed overlap gives a local minimum of the energy, and so a fixed point, near each pattern.
  - Zhou's result therefore gives fixed points near every ±1 pattern for all α < 0.1382. That covers T3, and it covers ±1 patterns too.
  - T3 stays a correct, elementary and self-contained proof for its model, but its constant is no longer a frontier.
- **"A ±1 bound above 0.113": known, if Zhou is correct.** The headline goal of this folder, a rigorous ±1 bound above 0.113 for fixed points near the patterns, is now claimed at 0.1382.
- **±1 universality as a technique: known in another form.** Zhou and Ali et al. prove universality for free energies and wells through Hamilton–Jacobi and Parisi-type formulas. Our proof of hypothesis (U) is different and elementary: capped Lindeberg plus a conditional Sudakov–Fernique comparison. It proves only what the barrier needs.
- **T2 (NLT barrier around every ±1 pattern for α < 0.12948): new as far as we found.**
  - The best rigorous NLT bound in print is 0.08 (Talagrand 1998).
  - Zhou's definition of a well allows m₂ = 1, the pattern itself, so his formula for the energy at fixed overlap does in principle determine the NLT capacity.
  - We found no statement or certificate of an NLT number in Zhou, Ali et al., or anywhere else.
  - T2 is the first rigorous proof of Stojnic's first-level NLT value for ±1 patterns, and it does not depend on any unrefereed preprint. It still needs its independent review before it can be called Proved.
- **T1: known in substance.** It is the Gaussian half of T2, and Zhou's machinery covers Gaussian rows.

## What is still open

1. **The exact NLT capacity for ±1 patterns.** It is the largest α for which some overlap m < 1 has a minimum energy above that of the pattern. Random duality theory predicts about 0.1298, so T2's 0.12948 is a lower bound and not the exact value.
   - Upper side: there is no rigorous upper bound on the NLT capacity below Zhou's AGS upper bound of 0.13820106.
   - Lower side: no proof goes past 0.12948.
   - A two-sided certified value is the natural next target. One route is certified evaluation of Zhou's zero-temperature formula at m₂ = 1, which would inherit that formula's unrefereed status. The other is a direct second-level lower bound.
2. **Independent checking of Zhou's 0.1382.** It is unrefereed and AI-written. Reproducing it would only verify an existing claim, so under this project's goal it is not a target.

## Retargeted work

- The headline of this folder changes from "beat 0.113" to "rigorous NLT barrier capacity for ±1 patterns ≥ 0.12948, up from 0.08".
- T3 is marked as superseded.
- Further work aims only at the open NLT capacity.
