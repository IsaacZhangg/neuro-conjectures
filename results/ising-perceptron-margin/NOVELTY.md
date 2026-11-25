# Novelty check (2026-11-25)

## Exact statement

For a margin κ > 0, let M_N(κ) be the largest M such that some J ∈ {−1,+1}^N satisfies
⟨g_μ, J⟩/√N ≥ κ for μ = 1, …, M, with i.i.d. standard Gaussian g_μ. The claim is that
M_N(κ)/N → α⋆(κ) in probability, where α⋆(κ) is the Krauth–Mézard replica-symmetric
threshold (the root of the RS free energy along its fixed-point branch).

Solving it at a given κ needs both directions there:
- upper bound: capacity ≤ α⋆(κ) (Huang's Theorem 3.6 once its four conditions hold);
- lower bound: capacity ≥ α⋆(κ) with high probability (a Ding–Sun type argument extended off κ = 0).

## How the search was done

Exa (web search and page fetch), WebSearch, and the Firecrawl research index (keyword search
from 2023 on, and the papers citing Ding–Sun 1809.07742, Huang 2404.18902 and Shmalo
2607.11010). Pages read in full or in part: Shmalo's GitHub README and the Zenodo pages of
all three versions of the margin preprint, Nakajima 2512.23195, Ding–Sun (journal version).

## What is known

- κ = 0 is settled: Ding–Sun (lower bound, conditional), Huang (upper bound, conditional),
  Shmalo arXiv 2607.11010 (all conditions certified in Arb), with the sharp threshold of
  Xu and Nakajima–Sun 2204.03469.
- Upper bound off zero margin: Shmalo, "The Ising perceptron at nonzero margin: certified
  bounds along the replica-symmetric curve", Zenodo 10.5281/zenodo.24181576 (v2.1, Nov 13,
  2026), code at github.com/yspennstate/ising-perceptron-margin. All four of Huang's
  conditions are certified at κ = −0.45, −0.05, 0.05, 0.0995 and 0.13, so capacity ≤ α⋆(κ)
  holds unconditionally there (for example α⋆(0.05) ∈ [0.781073068, 0.781074776]). Three of
  the four conditions (fixed point, AMP stability, local concavity) are certified on the
  strip κ ∈ [−0.65, 0.19]. The global two-variable condition (Huang's Condition 1.3) is
  certified only at the five margins.
- Uniqueness of the RS saddle point for every κ ≥ 0, analytically: Nakajima, arXiv 2512.23195.
  This covers the uniqueness part of Huang's Condition 3.1 for κ ≥ 0.
- Lower bound off zero margin: Shmalo's v1 (Zenodo 23329766) says it is "open at every
  nonzero margin"; it shows that unconditioned second moments fail at every margin, gives
  lower bounds for κ < 0 by embedding a symmetric perceptron (not stated as matching
  α⋆(κ)), and checks the general-margin Ding–Sun condition in floating point only at the
  five margins. Version 2.1 repeats: "Whether α⋆(κ) is also a lower bound at κ ≠ 0 is open."
- Non-rigorous or conditional: Stojnic, arXiv 2312.00073, recovers α⋆(κ) at every κ within
  his fully lifted random duality frame, assuming its stationarized second-level form.
- No paper or preprint found that proves the lower bound at any κ ≠ 0, or certifies all
  four upper-bound conditions at any κ > 0.13.

## What is open

1. The lower bound capacity ≥ α⋆(κ) at every κ ≠ 0. At κ ∈ {0.05, 0.0995, 0.13} a
   proof of it would settle the capacity outright, since the upper bound is already known.
2. The upper bound at κ > 0.13 (all four conditions), and three of the conditions past
   κ = 0.19.

## Verdict

Partly known. The κ = 0 case and the upper bound at κ ∈ {0.05, 0.0995, 0.13} are already
proved (Shmalo 2026). The lower bound at every κ ≠ 0 is open, and so is the upper bound
past κ = 0.13.

## What this means for the work here

- Earlier certifications in this folder of Huang's Conditions 3.1, 3.3, 3.4 and Claim 2.6 at
  κ = 0.1 are inside Shmalo's certified strip, so they are re-verification and are not
  claimed as new. The uniqueness part of Condition 3.1 for κ ≥ 0 is Nakajima's theorem.
- New relative to the literature found: Conditions 3.3, 3.4, the contraction and
  free-energy crossing of 3.1, and Claim 2.6 at κ = 0.25, 0.5, 0.75 and 1.0 (past κ = 0.19).
  These still leave Condition 1.3, so they are not yet an upper bound there.
- Retargeted work, in order:
  1. Lower bound at κ = 0.05 (then 0.13): write the κ = 0 steps of Ding–Sun for κ > 0 and
     certify the general-margin Ding–Sun condition in Arb. With Shmalo's upper bound this
     would determine the capacity at those margins. This is the main target.
  2. Upper bound at κ = 0.25 and 0.5: certify Condition 1.3, which with the conditions
     above gives capacity ≤ α⋆(κ) at margins nobody has reached.
