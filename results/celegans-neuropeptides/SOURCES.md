# Sources

Checked 2026-11-25. The literature screen was run by an AI agent; items it
could not open are listed under "Not verified".

## Primary

- Randi F, Sharma AK, Dvali S, Leifer AM. Neural signal propagation atlas of
  *Caenorhabditis elegans*. Nature 623, 406–414 (2023).
  doi:10.1038/s41586-023-06683-4 (preprint arXiv:2208.04790). Functional
  connection threshold q < 0.05 (Methods). Anatomy-constrained simulation
  predicts the atlas poorly in wild type and better in unc-31 animals
  (Fig. 3d; Extended Data Fig. 7). A screen for pairs connected in wild type
  but confidently not in unc-31 found 53 pairs (Extended Data Fig. 9); "the
  majority" of them express a matching peptide–receptor combination
  (Supplementary Table 1). That check is per pair, with no null or enrichment
  statistic.
- Ripoll-Sánchez L, Watteyne J, Sun H, et al. The neuropeptidergic connectome
  of *C. elegans*. Neuron 111, 3570–3589 (2023).
  doi:10.1016/j.neuron.2023.09.043. Short-, mid- and long-range models of the
  neuropeptide–GPCR network. We found no comparison with the signal-propagation
  atlas in the indexed full text.

## Later work that retests, refines or disputes it

- Sale C, Richardson (2026). Connectome structure constrains neural activity
  across scales. bioRxiv, doi:10.64898/2026.11.03.762695 (posted 2026-11-04).
  A Katz-type path predictor on the Cook 2019 wired connectome gives AUC 0.637
  over 28,989 atlas pairs and 0.603 (95% CI 0.566–0.640) among pairs with no
  direct synapse; a degree-preserving null drops it to about 0.517. Observation
  count alone gives AUC 0.709. Discussion: "the available monoamine and
  neuropeptide networks did not explain the measured responses". The details
  are in a supplement we could not read.
- Creamer MS, Leifer AM, Pillow JW (2024). Bridging the gap between the
  connectome and whole-brain activity in *C. elegans*. bioRxiv,
  doi:10.1101/2024.09.22.614271. A fitted connectome-constrained model;
  adding connections beyond the connectome did not significantly help on 56
  purely extrasynaptic pairs (p = 0.1, Fig. S4, per a summary of the PDF). No
  peptide network used.
- Moutuou & Benali, arXiv:2408.14221: adding single RID → URXL and RID → ADLR
  edges makes those responses significant in their model; the full peptide
  network is left for later work.
- Not about prediction with a peptide layer: Dvali, Seguin, Betzel, Leifer,
  arXiv:2412.14498; Sunil, Benali, Moutuou, arXiv:2604.02057; Reyner-Fuentes
  et al., arXiv:2610.37587; Moutuou & Benali, doi:10.64898/2026.10.29.762790;
  tool papers CeDNe (doi:10.1101/2025.11.03.683805) and the OpenWorm
  Connectome Toolbox (doi:10.64898/2026.11.04.764366).
- Cook SJ, Kalinski CA, Hobert O (2023). Neuronal contact predicts
  connectivity in the C. elegans brain. Curr Biol,
  doi:10.1016/j.cub.2023.04.071. Contact predicts synapses, not responses.
- Kinger, Bertram, Dyballa, Yemini, Zucker (2026), arXiv:2605.02852. Edges
  inferred from spontaneous activity, checked against synaptic, monoamine and
  GABA connectomes; no peptide layer with a proximity control.

## Datasets

- Signal-propagation atlas and default anatomical connectome (White et al.
  1986; Witvliet et al. 2021 adults), as bundled in the PyPI package
  wormneuroatlas 0.0.7.3 (wheel SHA-256 in fetch_data.py). Package code is
  GPL-3; we read only its data files and import none of its code. Data are not
  redistributed here.
- Ripoll-Sánchez et al. 2023 neuropeptide connectome CSVs
  (`01022024_neuropeptide_connectome_*_range_model.csv`), as bundled in the
  PyPI package cect 0.3.5 (MIT; wheel SHA-256 in fetch_data.py). Not
  redistributed here.
- Downloaded 2026-11-25 through PyPI. Zenodo, figshare, WormAtlas, OpenWorm
  and journal sites were unreachable from the cloud environment.

## Searches run

- Firecrawl research index: inspect and read for both primary DOIs, the Sale
  & Richardson DOI and the arXiv papers above; citers of both primary DOIs;
  six free-text searches (peptide connectome + signal propagation + prediction,
  extrasynaptic + AUROC, unc-31 + peptidergic edges, network communication
  models + stimulation, GNN or linear models + Randi atlas).
- Web search: seven queries on the same themes, plus the Sale & Richardson and
  Creamer titles.
- Read: PMC10632145 (Randi), ar5iv 2208.04790, Creamer PDF.

- Novelty audit, 2026-11-25: citers of both primary DOIs again; free-text
  searches since 2023 on peptide connectome versus function with a distance
  control, soma distance versus the atlas, and extrasynaptic connectome versus
  calcium imaging; re-read of the Sale & Richardson full text.

## Not verified

- Ripoll-Sánchez PhD thesis, doi:10.17863/CAM.115553 (PDF returned 503).
- Venkatesh et al. 2025 review, PMID 40135746 (no full text).

- Sale & Richardson supplement (fetch refused) and their code repository (not
  reachable from this session): which peptide network and statistic stand
  behind their negative sentence.
- Ripoll-Sánchez reference list (PMC served a browser check), so a passing
  citation of Randi cannot be excluded; which range model they treat as main.
- Exact numbers behind Randi Fig. 3d, Extended Data Fig. 7b and Supplementary
  Table 1.
