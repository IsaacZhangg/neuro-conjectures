# Novelty check

Checked 2026-11-25. Exa was unavailable (connection failure); searches used
Firecrawl paper search, Firecrawl's citation graph and web search instead.

## Verdict

**New.** We found no earlier rigorous proof that the all-pattern capacity of
dense associative memory with iid Gaussian patterns equals the union bound,
nor of the exact value (log 2)/4 at large λ. The primary source states this
as open. The lower bound (the union bound itself) is known and is not claimed
as new here.

## What was already known

- **Union-bound lower bound, all patterns, Gaussian.** Lucibello and Mézard,
  PRL 132, 077301 (2024), arXiv:2304.14964, eqs. (9)–(10) and App. F.2,
  eq. (86). Derived with the random energy model, physics-level. They write:
  "It remains an open question whether we have α_c^lb = α_c or not." Our
  Section 1 (lower bound) re-proves this bound rigorously; the value is
  theirs.
- **Typical-pattern threshold (one random pattern), Gaussian.** Lucibello and
  Mézard (2024); re-derived by Sumedha and Singh, arXiv:2607.19195 (2026),
  α_c = λ(1 − λ/2) for λ ≤ 1. Different question from ours (one pattern,
  not all).
- **Rigorous lower bounds for other pattern laws.** Demircigil, Heusel, Löwe,
  Upgang and Vermet, J. Stat. Phys. 168 (2017), arXiv:1702.01929: binary
  patterns, exponential capacity lower bound. Ramsauer et al., ICLR 2021,
  arXiv:2008.02217: patterns on a sphere, exponential capacity lower bound.
  Neither gives a matching upper bound or treats iid Gaussian patterns exactly.
- **Tight capacity for designed patterns.** Hu, Wu and Liu, arXiv:2410.23126
  (2024): optimal capacity when patterns form spherical codes. Not iid
  Gaussian.
- **Other related work, different question.** Petrova et al.,
  arXiv:2604.07401 and arXiv:2603.13350 (spherical patterns, typical
  retrieval); Ota and Taki, arXiv:2609.10976 (phases, typical Gaussian
  patterns); Agliari et al., arXiv:2607.29554 (binary, multilayer, annealed,
  one-sided); Barnfield et al., arXiv:2605.05189 (sharp top-1 threshold for
  linear associative memory, a different model); Albanese et al.,
  arXiv:2604.02789 and arXiv:2509.06905; Clark, arXiv:2506.05303; Achilli et
  al., arXiv:2502.09578.

## What is new here

- A matching **upper bound**: for α above the union bound, with probability
  tending to 1 some pattern fails the static retrieval criterion (PROOF.md
  Sections 2–4, second-moment arguments). Hence α_c(λ) = α_c^lb(λ) for every
  λ > 0, including λ = ∞, under LM's static criterion.
- The **exact large-λ value** α_c = (log 2)/4 for λ ≥ λ_x = 0.700917829…,
  with λ_x certified by interval arithmetic (check.py).
- A **sharp transition**: below α_c all patterns are retrieved with
  probability tending to 1, above it with probability tending to 0.

## Not covered

- The dynamical definition of retrieval (convergence of the update rule from
  a perturbed start). Only LM's static signal-dominance criterion is proved.
- Possible classical precursors of the λ = ∞ geometric event (every Gaussian
  vector has larger squared norm than its inner product with every other).
  We searched spherical-cap and random-geometry work (e.g. arXiv:1511.06198,
  arXiv:1612.09045, arXiv:2610.20299) and found nothing stating this
  threshold, but cannot rule out an older source in that literature.

## Searches run

- Firecrawl paper search, 2023 onward: dense associative memory capacity,
  Gaussian patterns, all patterns, union bound, rigorous, second moment,
  exponential Hopfield.
- Citation graph: papers citing arXiv:2304.14964, 30 screened by title and
  abstract, relevant ones read.
- Rigorous-capacity search: Demircigil et al., Löwe and coauthors, Ramsauer
  et al. and their citers.
- Web search for "all patterns", "union bound", "log 2", "Lucibello Mézard
  open question".
- Geometric search on the λ = ∞ event (maximal inner products of Gaussian
  vectors, spherical caps).
