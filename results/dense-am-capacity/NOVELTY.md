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
- Nothing beyond LM's static criterion. The λ = ∞ case is checked against
  the random-packing literature in the next section.

## The λ = ∞ case and random spherical codes

At λ = ∞ the criterion is that ‖ξ^μ‖² > ξ^μ·ξ^ν for every pair μ ≠ ν,
i.e. no pattern has a larger inner product with another pattern than with
itself. We checked whether the (log 2)/4 threshold already follows from the
extreme-angle results for exponentially many random points:

- Cai, Fan and Jiang, Distributions of angles in random packing on spheres,
  JMLR 14 (2013), arXiv:1306.0256, Theorem 6: for n uniform points on the
  sphere in R^p with (log n)/p → β ∈ (0, ∞), the minimum pairwise angle
  converges to arccos √(1 − e^{−4β}).
- Cai and Jiang, Phase transition in limiting distributions of coherence of
  high-dimensional random matrices, J. Multivariate Anal. 107 (2012),
  arXiv:1102.2926: the same regime for the largest Pearson correlation
  between columns.

Verdict for λ = ∞: **does not follow.** Both results are about normalized
quantities (angles, correlations), which do not see the pattern norms. For
points of equal norm the λ = ∞ criterion holds automatically, since
cos θ < 1 for distinct points, so the Gaussian threshold comes entirely from
norm fluctuations acting together with the angles: a short pattern
(‖ξ‖²/N ≈ 1/2) loses to one competitor. Neither paper treats that event,
and the angle limit alone does not determine it.

There is a numerical coincidence worth recording. Cai, Fan and Jiang's
minimum angle equals π/4 exactly at β = (log 2)/4, our threshold. The reason
is that both pair events have the same large-deviation rate per dimension,
(log 2)/2: P(cos θ > 1/√2) ≈ 2^{−N/2}, and P(ξ·η > ‖ξ‖²) ≈ e^{−N min_s
[I(s) + s/2]} = 2^{−N/2} (minimum at s = 1/2). The events and their
second-moment structure differ, so neither result implies the other; our
second-moment upper bound is in the same spirit as their exponential-case
argument. We also screened 40 papers citing Cai, Fan and Jiang (including
Kai Zhang, arXiv:1306.0623 and arXiv:1511.06198, on maximal inner products
of normalized vectors) and found none treating the unnormalized event.

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
  vectors, spherical caps, maximum inner product search).
- Cai, Fan and Jiang (2013) and Cai and Jiang (2012) read; papers citing
  Cai, Fan and Jiang screened by title and abstract.
