# Methods

## Model

- N binary neurons s_i ∈ {−1, +1}.
- P = round(αN) patterns, each entry ±1 with probability 1/2, independent.
- Couplings W = ξᵀξ / N with zero diagonal (no self-coupling).
- Local field h = W s.

## Dynamics

Zero-temperature asynchronous updates. At each step, one neuron is chosen
uniformly among the neurons whose sign disagrees with their field and is
flipped; the field is updated incrementally. Each flip lowers the energy, so the
dynamics stops at a fixed point. Choosing only among unstable neurons gives the
same set of reachable fixed points as random sequential sweeps, with less
computation.

## Design

- Sizes N = 500, 1000, 2000, 4000.
- Loads α = 0.10 to 0.24 in steps of 0.01.
- 12 independent networks per (N, α), 10 starting patterns per network,
  120 samples per point.
- Single NumPy generator, seed 20261125.

## Analysis

As defined in CLAIM.md. The script exits 0 if the decision rule passes.

## Robustness checks to add before review

- Double the samples per point and confirm the interpolated crossings move by
  less than 0.005.
- Repeat with synchronous updates and report the difference.
- Try a 1/N correction instead of N^-1/2 and report the spread of intercepts.

## Deviations from plan

- A first quick run used α up to 0.20 and put the N = 250 crossing at 0.18,
  close to the edge of that grid. Before the full run, the grid was widened to
  0.24 and the full-run sizes were set to 500 to 4000. No change was made after
  the full run.
