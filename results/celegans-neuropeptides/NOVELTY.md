# Novelty check

Checked 2026-11-25 by an AI agent. Exa was unavailable (connection failure);
searches used the Firecrawl paper index, its citation graph (papers citing
Randi et al. 2023 and Ripoll-Sánchez et al. 2023) and Firecrawl web search,
covering 2023 to now.

## Verdict

**Partly known.** The wired-connectome baseline and the coarse "peptides do not
explain the atlas" answer are already in print. The pre-registered,
peptide-specific test, the gain on pairs with no wired path, and the
proximity control that explains that gain are not, as far as we found.

## What was already known

- **Peptides line up with some wireless responses, pair by pair.** Randi et
  al., Nature 623, 406–414 (2023), doi:10.1038/s41586-023-06683-4: of 53
  pairs connected in wild type but not in unc-31, "the majority" express a
  matching peptide–receptor pair (Extended Data Fig. 9, Supplementary
  Table 1). No null, no network-wide statistic.
- **Path score on wiring predicts the atlas; degree-preserving null.** Sale
  and Richardson, bioRxiv doi:10.64898/2026.11.03.762695 (posted 2026-11-04):
  a parameter-free path operator on the Cook 2019 connectome gives AUC 0.637
  over all observed pairs and 0.603 among unsynapsed pairs, with most of the
  gain from two-step paths, scored against degree-preserving surrogates. Our
  wired baseline, S(A) = A + A·A with a rewiring null, is the same idea,
  frozen independently; **that part is not new.**
- **Peptide networks do not explain the responses, in one sentence.** Same
  paper, Discussion: "the available monoamine and neuropeptide networks did
  not explain the measured responses". The network, statistic and any spatial
  control sit in a supplement we could not read. Our "no for the network as a
  whole" agrees with this and is **not new as a headline**.
- **Fitted models need no extra edges.** Creamer, Leifer and Pillow (2024),
  doi:10.1101/2024.09.22.614271: adding connections beyond the connectome did
  not significantly help on 56 purely extrasynaptic pairs. No peptide network.
- **Single peptide edges in a model.** Moutuou and Benali, arXiv:2408.14221:
  adding RID → URXL and RID → ADLR makes those two responses significant.
- **Contact predicts synapses.** Cook, Kalinski and Hobert, Curr. Biol. (2023),
  doi:10.1016/j.cub.2023.04.071: neurite contact predicts synaptic
  connectivity. About wiring, not about measured responses or peptides.
- **Lag-specific edges from activity.** Kinger, Bertram, Dyballa, Yemini and
  Zucker, arXiv:2605.02852 (2026): edges inferred from spontaneous activity,
  scored against synaptic, monoamine and GABA connectomes, with the Randi
  atlas as a benchmark. Not a peptide-network test with a proximity control.

## What is new here

1. A pre-registered test (CLAIM and METHODS frozen before data) of whether
   adding the Ripoll-Sánchez peptide connectome to wiring improves AUROC on
   the Randi atlas, with a cluster bootstrap over stimulated neurons and a
   peptide-specific degree-preserving null, split by pairs with and without a
   wired path.
2. The positive part: a small gain on pairs with no wired path, short-range
   model only.
3. The proximity follow-up (frozen before it ran): soma nearness alone gives
   AUROC 0.72 on no-wired-path pairs; within matched distance and ganglion
   classes the pair-specific peptide gain is +0.003 (CI −0.042 to +0.051), so
   **the gain is mostly explained by physical proximity.** We found no earlier
   test of the peptide connectome against the atlas with a proximity control.

## Not verified

- Sale and Richardson supplement: which peptide network, statistic and
  controls stand behind their negative sentence. If it contains a
  proximity-matched test, item 3 above would be partly known too.
- Ripoll-Sánchez PhD thesis, Cambridge 2025, doi:10.17863/CAM.115553 (PDF
  returned 503). Its abstract says it looked at how neuropeptides interact
  with synaptically wired circuits; we could not see whether it compares with
  the atlas.
- Venkatesh et al. 2025 review (PMID 40135746): no full text reached.
