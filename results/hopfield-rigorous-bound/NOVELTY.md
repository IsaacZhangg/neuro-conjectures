# Novelty audit

Literature search done on 2026-11-25. Tools: Exa web search, the Firecrawl research index (arXiv, PubMed, bioRxiv) from 2020 to now, and papers citing Stojnic (2024). arXiv's own site turned away direct requests, so arXiv papers were read through these indexes. Google Scholar was not reachable from here.

**Verdict: partly known.** Most of what this folder set out to do was done by a preprint posted on 2026-11-05. That preprint has not been refereed. One part of this folder still looks new: the rigorous bound for the firm energy barrier around each ±1 pattern.

**The barrier bound and Zhou's capacity measure different things.** Zhou's 0.1382 is the AGS capacity: a load has an energy well somewhere near each pattern. Our 0.12948 is the NLT barrier capacity: every configuration at Hamming radius ⌈δN⌉ has higher energy than the pattern itself, so zero-temperature dynamics started at the pattern stay nearby. The barrier capacity is at most the AGS capacity, so the two numbers are not competing values for one quantity, and 0.12948 is not a weaker version of 0.1382.

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

## Second check (2026-11-25, after reading Zhou in full)

- **Read:** about two thirds of the text of hexagon 2611.00105v1, word for word.

### Zhou's notion of capacity

- **The notion is an AGS well in the energy at fixed overlap.** At T = 0, a load α has a well if e(m₂) > e(m₁) for some m_s ≤ m₁ < m₂ ≤ 1, where e is the limiting maximum of −H/N at fixed overlap.
- **m₂ = 1 is allowed but not required.** We found no value or certificate for m₂ = 1, the energy barrier around the pattern itself.
- **NLT is never mentioned.** The text names neither NLT nor Talagrand's 0.08 nor Stojnic's 0.1295 or 0.1298.
- **Both pattern models are covered.** The paper covers ±1 patterns and Gaussian rows, and proves they have the same wells.
- **It concerns minima of fixed-overlap slices, not dynamics.** It holds for all patterns at once (Theorem 1.2(e)). It says nothing about the dynamics. Its §1.6 sets finite-N retrieval dynamics apart from the wells of the limiting profile.

### Inferences (ours, not stated by Zhou)

- **Fixed points:** a well of this kind gives a single-flip local minimum near each pattern. One flip moves the overlap by 2/N, so the minimizer over the overlap band above m₁ is stable to single flips.
- **Single-flip, every-pattern version: covered by Zhou, if correct.** "Every pattern at once has a single-flip-stable state nearby" holds up to 0.1382 by this inference. It is not open.

### Soundness

- **Unrefereed dependence.** The paper rests on Ali et al., arXiv:2611.03306, which is unrefereed.
- **Computer-assisted parts.** Theorem 1.2 (c), (d) and (f) are computer-assisted.
- **Not re-verified.** We did not re-verify the paper, because pure re-verification is out of this project's scope. We found no error in what we read, but that is not a check.

### Superseded or known, if Zhou is correct

- **The original goal.** A rigorous ±1 bound above 0.113 for a stable state near each pattern.
- **Theorem 3.** A fixed point for α < 0.134716, Gaussian interference.
- **The fixed-point reading of Theorems 1 and 2.** As fixed-point results they sit below 0.1382.

### Still open and not covered by Zhou: retrieval from the pattern

- **The energy barrier around the pattern itself (NLT).** This is the barrier H(σ) > H(ξ) on the whole sphere at radius ⌈δN⌉.
- **What it guarantees.**
  - This barrier is what makes zero-temperature single-flip dynamics, started at the pattern, stay within ⌈δN⌉ of it. The energy never goes up, so the dynamics cannot cross the sphere.
  - An AGS well alone does not give this. The well near capacity can sit above the pattern's own energy, so the dynamics started at the pattern may leave.
- **Where the bounds stand.**
  - The best rigorous ±1 value in print is 0.08 (Talagrand 1998).
  - This folder's Theorem 2 proves 0.12948, with its review pending.
  - Random duality theory predicts the true value is about 0.1298.
  - No rigorous upper bound exists below 0.1382.
- **Finite-N rates.**
  - Our bound holds with probability 1 − 5P·e^{−cN}, with c not explicit.
  - Making c explicit would be an extension, but only a minor one.

### Decision

- **Keep:** the claim of novelty only for the dynamical retrieval version, the ±1 energy barrier at 0.12948.
- **Finish:** its independent review.
- **Don't spend a new slot on:** pushing it from 0.12948 to the predicted 0.1298. The gain is small, and the clean route goes through Zhou's unrefereed formula.
