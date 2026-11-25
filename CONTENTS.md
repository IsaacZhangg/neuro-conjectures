# Catalogue

Every result with its exact scope. The README table is the short form.

## Model claims

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Hopfield network storage capacity](results/hopfield-capacity/) | Hebbian Hopfield network, ±1 neurons, unbiased random patterns, zero-temperature asynchronous dynamics started at a stored pattern; N = 500 to 4000; capacity defined as the load where the probability of ending with overlap ≥ 0.9 falls below one half; extrapolated to large N with an assumed N^-1/2 correction. | Simulation, fixed seed, recorded output. Replication of Amit, Gutfreund and Sompolinsky (1985). Not yet independently reviewed. |

## Theorems

None yet.

## Empirical claims

| Result | Exact scope tested | Evidence |
| :-- | :-- | :-- |
| [Cuntz 2/3 wiring law on EM connectomes](results/cuntz-wiring-law/) | Branch-point form L ∝ V^(1/3)·bp^(2/3) of Cuntz, Mathy and Häusser (2012); slope of log(L/V^(1/3)) on log(bp) with convex-hull V and 1 µm twig pruning; FlyWire 783 whole-neuron skeletons (402 cell types, 15 neurons each, one fly) and MICrONS minnie65 v661 dendritic skeletons (1,309 E/I cells of the V1 column, one mouse). FlyWire Refuted (β = 0.537 [0.532, 0.541]); MICrONS Contested (β = 0.528 [0.514, 0.541], but 0.613 with tips as n). | Pre-registered reanalysis, plan committed before the data were read, recorded output. Not yet independently reviewed by a person. |

Screened targets are in [docs/candidates.md](docs/candidates.md).
