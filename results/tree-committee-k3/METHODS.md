# Methods

## Model

Tree committee machine, K = 3, n = 3N inputs in three blocks of N, output
sign(Σ_k sign⟨w_k, ξ_k⟩). Patterns have i.i.d. N(0,1) entries; m = ⌊αn⌋, so
each branch sees load β = m/N = 3α. After multiplying each pattern by its label,
storing everything means: every pattern has at least two branches with positive
field.

## How the bound is built

1. **Search for a route (exploratory).** A Cover-only construction gives 4/3.
   Replacing the first two branches by convex programs and computing their
   error fractions with the convex Gaussian min-max theorem gives more.
2. **Choice of loss.** The first attempt used plain hinge loss and assumed the
   weight-norm constraint binds; simulation (error 0.273 against a predicted
   0.193) showed the assumption was wrong. The fix is the binding condition
   P = Σ L E[(G + D)D] > 0, and a linear term −ρh in the loss, which keeps P
   positive while lowering the error fraction. With it, theory and simulation
   agree (predicted 0.2266, simulated 0.2228 at κ = 0.02).
3. **Constants.** κ = 1/1000, ρ₁ = 59/250, ρ₂ = 7/25, padding a = 975253/10⁶,
   chosen numerically so that the stage-3 load stays below 2 with P bounded away
   from 0 at both stages. ρ₂ = 0.279 failed the certified P > 0 test, so 0.28 is
   used. The best α found this way is about 1.5517; the theorem uses 1.55.
4. **Proof.** PROOF.md: exact penalty for the hard constraints of stage 2
   (multipliers bounded through a Slater point), lifted convex Gaussian min-max
   theorem, convexity of the auxiliary problem in the norm, strong convexity of
   the Lagrangian for the deviation step, and a Lipschitz sandwich for the error
   count.

## Certified numerics

check.py uses python-flint (Arb ball arithmetic, 200 bits). All Gaussian
expectations are closed forms in Φ and φ. For each stage it brackets the root of
c(x) = 1 to width 10⁻⁶ by bisection on certified signs, encloses P on the whole
bracket, and checks six inequalities: P > 0 (both stages), the stage-1 error
load below the padding a, a < 1, the Slater condition 2κ√a < 1 − √a, and the
stage-3 load a + (β − a)Φ(κ − x₂) < 2. Recorded result: stage-3 load
1.99950566571333 ± 3.2·10⁻¹⁵.

## Sanity check by simulation

sim_pipeline.py runs all three stages at finite N with cvxpy (Clarabel) and
κ = 0.02. Seed 0 stored every pattern at α = 1.4 (N = 250; stage-3 load
1.74 N) and α = 1.5 (N = 400; stage-3 load 1.83 N), rerun 2026-11-25 with
`python sim_pipeline.py 250 1.4 0` and `python sim_pipeline.py 400 1.5 0`.
Error fractions were close to the limit predictions (stage 1: 0.179 and 0.198
against 0.198 and 0.212; stage 2: 0.287 and 0.259 against 0.250 and 0.276);
the stage-2 excess at N = 250 is larger than sampling noise alone, which
points to finite-N corrections; it is gone at N = 400. This is not part of the proof; it checks that the limit formulas
describe the optimisers.

## Robustness

- The margin κ only needs to be positive (it keeps the solution away from
  w = 0 at the boundary of the binding condition); 1/1000 costs little.
- Leaving some patterns untried in stage 1 or 2 (they then pass with probability
  1/2, independently) did not improve the bound.
- The non-convex ideal of the same three-stage scheme (Gardner–Derrida margin
  laws per stage, not rigorous) gives about 2.41, which marks how much the
  convexity restriction costs.

## Deviations from plan

The first numeric estimate (1.616) rested on the wrong binding assumption in
step 2 and was discarded. No other deviations.
