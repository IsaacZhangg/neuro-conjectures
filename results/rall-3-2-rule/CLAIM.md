# Claim

## Statement tested

At a dendritic bifurcation with parent diameter d₀ and daughter diameters d₁
and d₂,

  d₀^(3/2) = d₁^(3/2) + d₂^(3/2),

that is, the branch exponent η in d₀^η = d₁^η + d₂^η equals 3/2.

This is our reading of the question "do real dendrites satisfy Rall's 3/2
power rule?". Rall stated the relation as a condition under which a branched
tree maps onto an equivalent cylinder (together with equal electrotonic
lengths to all terminals), not as an empirical law. The empirical claim tested
here is the one later work tested against measured trees (Rall et al. 1992 as
summarised by Rall 2009; Desmond and Levy 1984; Larkman et al. 1992): that the
condition approximately holds at the bifurcations of real dendrites.

## Where it comes from

- Rall (1959), Exp. Neurol. 1:491–527, where the d^(3/2) dependence of branch
  input conductance leads to the constraint. Statement location in the full
  text not confirmed (paywalled; see SOURCES.md).
- Rall (2009), Scholarpedia 4(4):1369, "Rall model", section "d3/2 Constraint
  on branch diameters": "the sum of their 3/2 power values is equal to the 3/2
  power of the parent branch diameter. When this constraint is satisfied at all
  branch points of a tree, this tree can be mapped to an equivalent cylinder
  (provided also that all terminal branches end at the same electrotonic
  distance from the trunk)."

## Operational definitions

- Dendrite: SWC samples of type 3 (basal dendrite) or 4 (apical dendrite) in
  the NeuroMorpho.org standardised (CNG) files.
- Bifurcation: a dendritic node with exactly two dendritic children and a
  dendritic parent (not the soma).
- d₀, d₁, d₂: mean diameter of the samples within 5 µm path distance of the
  bifurcation node on the parent side (including the node) and on each
  daughter side (excluding the node), stopping at the next branch point or
  the soma. If no sample falls in a window, the nearest sample is used.
- Per-cell exponent η̂: the value in [0.5, 6] that minimises
  Σ_b [ (1/η)·ln(d₁^η + d₂^η) − ln d₀ ]² over the cell's bifurcations b.
- Group: NeuroMorpho archive (one lab) × species × cell-type label.
- Group estimate: median of per-cell η̂, with a 95% percentile bootstrap
  interval over cells (2000 resamples).
- Pooled estimate: median of group medians, with a 95% percentile bootstrap
  interval over archives (2000 resamples, all of an archive's groups resampled
  together).
- Margin: η ∈ [1.35, 1.65]. At a symmetric bifurcation this is a daughter to
  parent diameter ratio of 0.60 to 0.66 against 0.63 for η = 3/2, about the
  diameter precision of light-microscope reconstructions. The margin is our
  choice, fixed here before any data was read.

## Decision rule (committed before any code reads NeuroMorpho data)

Each qualifying group is labelled:

- consistent if its 95% interval lies inside [1.35, 1.65];
- deviating if its 95% interval does not overlap [1.35, 1.65];
- inconclusive otherwise.

Overall status, from the primary analysis:

- **Supported** if the pooled 95% interval lies inside [1.35, 1.65] and at
  least two thirds of qualifying groups are consistent.
- **Refuted** if the pooled 95% interval does not overlap [1.35, 1.65] and at
  least two thirds of qualifying groups are deviating.
- **Partial** otherwise, if at least one group is consistent; the README then
  names the groups where the rule holds and where it fails.
- If fewer than 5 groups qualify, or no group is consistent and the Refuted
  condition is not met, the result stays **Open** and is reported as
  inconclusive.
- **Contested** if any robustness analysis in METHODS.md gives a different
  status from the primary analysis.

`check.py` exits 0 only when the primary status is Supported and no robustness
analysis changes it, and 1 otherwise. The status word is written to
`check-output.json`.
