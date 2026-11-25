# Sources

Checked 2026-11-25.

## Primary

- Jeong H, Taylor A, Floeder JR, Lohmann M, Mihalas S, Wu B, Zhou M, Burke DA,
  Namboodiri VMK. Mesolimbic dopamine release conveys causal associations.
  Science 378, eabq6740 (2022). doi:10.1126/science.abq6740. Tests 1 to 11:
  Figs 3 to 6, S8 to S12, Supplementary Notes 6, 9, 10. Analysis and simulation
  code: github.com/namboodirilab/ANCCR (`analysis/`, `simulation/`).

## Later work that retests, refines or disputes it

- Loewinger G et al. 2025, eLife, doi:10.7554/eLife.95802. Reanalysis of the
  same data with functional mixed models: the Test 1 rise holds across sessions
  but responses fall within a session; Test 2 holds; Test 4 is challenged (the
  late part of the cue response is lower after the delay extension).
- Gershman SJ et al. 2024, argues verbally that TD with a rich state space can
  produce Test 1; no simulation of all tests.
- Qian L et al. 2025, TD with inter-trial-interval states explains a
  contingency-degradation effect like Test 7, in their own task.
- Amo R et al. 2022, Nature Neuroscience: backpropagation of dopamine during
  learning in a different task (in tension with Test 9).
- Cone I, Clopath C, Shouval HZ 2024, Nature Communications: discusses Tests 9
  and 10 and learning without backpropagation.
- Burke DA et al. 2026 and Jeong H et al. 2026 (memory erasure): further tests
  that the authors report TD fails.

None of these runs a single TD model through all eleven tests on the animals'
event logs, and none proves that a test is unreachable for a defined TD family.

## Datasets

- DANDI 000351 version 0.260828.1905, "Mesolimbic dopamine release conveys
  causal associations", CC-BY-4.0, files read from the DANDI S3 bucket
  2026-11-25.

## Searches run

- Exa was not reachable from this session; Firecrawl research search and web
  search were used for: "Jeong 2022 dopamine causal associations TD model",
  "ANCCR temporal difference reanalysis", "Jeong 2022 reanalysis dopamine",
  "TD learning reproduce ANCCR tests", "backpropagation dopamine Jeong", and
  citation searches from the primary paper, 2023 to 2026.
