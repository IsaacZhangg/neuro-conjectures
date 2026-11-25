# Claim

## Statement tested

Target #65 from the conjecture finder's list, frozen here word for word:

> No temporal-difference (TD) reward prediction error model with one fixed state
> representation and one parameter set per animal can reproduce the full set of
> nucleus accumbens dopamine results in Jeong et al. 2022 that the retrospective
> causal-learning model ANCCR reproduces.

Counts as solved if (frozen word for word):

> Either a single TD model (fixed state representation, one parameter set per
> animal) is exhibited that passes every experiment in the public dataset, which
> disproves the statement, or it is shown by proof or by an exhaustive argument
> over a defined family of TD models that a named experiment cannot be
> reproduced.

This is an open question, not a re-verification. Jeong et al. 2022 showed that
the specific TD models they simulated (TD(0)/TD(λ) with complete serial compound
or microstimulus states, trial-based resets) fail their tests. As of the status
check on 2026-11-25 (see SOURCES.md) nobody has run one TD model through all
eleven tests on the animals' own event logs, and nobody has proved that a test
is out of reach for a defined TD family.

## Where it comes from

Jeong H, Taylor A, Floeder JR, Lohmann M, Mihalas S, Wu B, Zhou M, Burke DA,
Namboodiri VMK. Mesolimbic dopamine release conveys causal associations.
Science 378, eabq6740 (2022). doi:10.1126/science.abq6740. Tests 1 to 11 are in
Figs 3 to 6, S8 to S12 and Supplementary Notes 6, 9 and 10. The paper says
(Test 1 reasoning) that a TD RPE to a reward "cannot be higher than that for an
unpredicted reward", and treats Tests 1 to 11 as predictions TD RPE gets wrong.

## Our reading (operational)

"TD model" means a model in the family TD-MS defined in METHODS.md: linear value
over a fixed feature map of the event history, TD(λ) learning, dopamine equal to
the TD error δ seen through a fixed photometry kernel. "One fixed state
representation" means the feature map and its hyperparameters do not change
across sessions or protocols. "One parameter set per animal" means each animal
gets one parameter vector for its whole history; we try first one vector shared
by every animal, which is stricter.

"Reproduce" means: the model is driven by each animal's real event log from
DANDI 000351, every session in order, and its simulated dopamine is analysed
with Jeong's own measures and statistics (METHODS.md, section "Tests"). A test
passes when the model reaches the paper's conclusion by the pass rule written for
that test. Licking is not modelled: every test that compares dopamine with
behaviour uses the animal's real licks.

## Decision rule (committed before any code reads the data)

- **Refuted** (the statement is disproved, which solves the target): one model
  in TD-MS, with one parameter vector per animal, passes all eleven tests by the
  primary rule, and still passes them under the robustness checks in
  METHODS.md (other noise seeds, the other missing-day rule, dt = 0.05 s).
- **Supported** (the target is solved the other way): a written argument shows
  that no model in a stated subfamily of TD-MS can pass a named test. Only a
  proof counts here. A parameter search that finds nothing is reported as
  "not found", never as Supported.
- **Contested**: a model passes all eleven tests by the primary rule but a
  robustness check flips at least one test.
- **Partial**: neither of the above. We report the best model, how many tests it
  passes, which ones fail, and any proof that covers part of the family.
