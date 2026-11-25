# Novelty audit for target #65

Audit date: 2026-11-25. Asked for by the project's coordinator, quoting Isaac:
"Make sure for all of the conjecture solvers channels that none of them are
proving things that have already been proved, and that we're trying to push new
boundaries and make new discoveries and prove unproved conjectures."

**Verdict: partly known.** Single tests have TD accounts in print, mostly on
other datasets. Nobody has fit one TD model to all eleven Jeong tests on
Jeong's own data, and nobody has proved any test out of reach for a defined TD
family. That joint question, which is the frozen claim, is open. The target
stays as frozen in CLAIM.md; no test or class is dropped.

## How we searched

- Exa web search on TD models of Jeong 2022, belief-state and successor
  TD for contingency degradation, and rebuttals of Jeong 2022.
- Firecrawl paper index (PubMed, bioRxiv, medRxiv, arXiv), two queries on TD
  models and ANCCR / Jeong 2022. Its citation-graph lookup for Jeong 2022
  returned nothing, so citing papers were found by search instead.
- Main text read for Qian et al. 2025, Loewinger et al., Zhang et al. 2025
  and the Uchida 2024 essay. Amo et al. 2022 is cited from prior knowledge and
  was not re-read here.
- Window: 2022 to 2026-11-25.

## What is known

| Source | What it shows | Data | Bears on |
|---|---|---|---|
| Qian, Burrell, Hennig, Matias, Murthy, Gershman, Uchida, *Nature Neuroscience* 2025, doi:10.1038/s41593-025-01915-4 (bioRxiv 10.1101/2024.02.05.578961) | TD with a belief state over the inter-trial interval (or a cue-context state) explains the drop in cue-evoked dopamine under contingency degradation. Their microstimulus TD did not (Extended Data Fig. 5). ANCCR failed their cued-reward control. | Their own odor task, not Jeong's | Test 7, for belief-state TD, in principle |
| Zhang, Costa, Langdon, Schoenbaum, *Trends in Cognitive Sciences* 29:434 (2025), doi:10.1016/j.tics.2025.02.001 | Microstimulus TD with separate basis sets for each cue occurrence gives a positive but weaker intermediate-cue response and a larger reward response after it, matching Jeong's trial-less ratios (Box 2). Simulation only, code on Zenodo. | Task layout copied from Jeong; not fit to Jeong's recordings | Test 8, for TD-MS with per-occurrence features |
| Loewinger et al., *eLife* (10.7554/eLife.95802) | Functional mixed-model reanalysis of Jeong's own photometry. The positive trend of reward response over random-reward sessions is a Simpson's paradox: within-session slopes are negative, between-session intercepts rise. The delay-change test shows lower dopamine 1 to 2.5 s after the cue on long-delay trials. The backpropagation result agrees with Jeong. | Jeong's data | Tests 1, 4, 9 (statistics, not a model) |
| Uchida, *The Transmitter* essay, 2024 | Argues Jeong tested only simple TD variants, and that TD with suitable parameters or representations explains part of the data (citing Qian, and Pan et al. 2005 for sequential conditioning). No new fit. | none | Framing |
| Cone, Clopath, Shouval, FLEX model (PMC10543312) | A plastic-time-representation model where CS1 and CS2 responses co-evolve, like Jeong's sequential result. Not a TD RPE model with fixed features. | Simulation | Test 10, for a non-TD model |
| Gershman et al., *Nature Neuroscience* 2024, "Explaining dopamine through prediction errors and beyond" | Review; lists Jeong 2022 among open challenges to TD. | none | Framing |
| Amo et al., *Nature Neuroscience* 2022 | Reports gradual backward shift of dopamine during learning, as TD predicts. | Their own task | Test 9 conflict across datasets, not on Jeong's |
| ANCCR side: Garr et al. 2024 (*Science Advances*); Burke et al. 2026 (*Nature Neuroscience*, doi:10.1038/s41593-026-02206-2); Jeong et al. 2026 bioRxiv (doi:10.64898/2026.01.17.700100) | Further ANCCR tests on new tasks. None fits TD to the 2022 data. | New datasets | Context |

## What is open

1. One TD model, one parameter set per animal, run on Jeong's event logs and
   scored on all eleven tests. Every paper above handles one test, on its own
   task or a copied layout. None checks whether the fixes are compatible.
2. A proof that some test is out of reach for a defined TD family. None exists
   in either direction.
3. Whether the fixes conflict. Our interim screen (2 noise seeds, Tests 1 to 9)
   already suggests they do:
   - TD-MS (883 points): Test 7 never passes; Test 8 passes rarely, which fits
     Zhang et al.
   - TD-BS (355 points): Test 7 passes, which fits Qian et al.; Test 1 never
     passes.
   No point in either class has passed more than six of nine. A
   representation that helps Test 7 seems to lose Test 1, and the reverse.
   That trade-off is not in print.

## What this changes in our plan

- Nothing is dropped. Test 7 for TD-BS and Test 8 for TD-MS each have a
  precedent, so passing either alone is not a finding. The finding is the
  joint pass, or a proof that it cannot happen.
- Loewinger's Test 1 reanalysis is noted, but the frozen pass rule for Test 1
  stays as written (AGENTS.md: never reword claims). It is reported beside the
  result.
- The impossibility work moves toward the trade-off in point 3: whether a
  TD-MS or TD-BS subfamily can pass Tests 1 and 7 together.
