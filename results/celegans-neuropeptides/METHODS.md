# Methods

Committed before any code reads the data. The git history is the record.
Until this file was committed, only file listings, file sizes, checksums and
the loader code of the two packages were inspected, to learn file formats and
index conventions. No data values were read. check.py was tested only on a
synthetic network before this commit.

## Data

All files come from two PyPI packages, pinned by SHA-256 in `fetch_data.py`.
No package code is imported; the files are read directly.

| Item | File | Source package |
| :-- | :-- | :-- |
| Signal-propagation atlas (wild type and unc-31) | `funatlas.h5` | wormneuroatlas 0.0.7.3 |
| Neuron order of the atlas | `neuron_ids.txt` (second column) | wormneuroatlas 0.0.7.3 |
| Wired connectome | `aconnectome_default.h5` (`chem`, `gap`) | wormneuroatlas 0.0.7.3 |
| Neuropeptide connectome | `01022024_neuropeptide_connectome_{short,mid,long}_range_model.csv` | cect 0.3.5 |

Conventions taken from the package loaders:

- Atlas and anatomy matrices are indexed `[post, pre]`: entry (i, j) is the
  effect of j on i. Groups `wt` and `unc31` each hold `dFF`, `q`, `q_eq` and
  `occ1` (number of trials).
- The wired connectome is the default the atlas package uses: chemical and
  gap-junction counts summed over White et al. 1986 and Witvliet et al. 2021
  adult datasets.
- Peptide CSV: first row and first column are neuron names; row = releasing
  (pre) neuron, column = receiving (post) neuron; entries count
  neuropeptide–GPCR pathways. Names are standardised as in cect (DA01 → DA1,
  and the same for DB, AS, VA, VB, VC, VD, DD). AWCL → AWCOFF and AWCR → AWCON,
  the convention the atlas package uses. Any peptide-network neuron that still
  has no atlas match is dropped and listed in the output; atlas neurons with no
  peptide-network row get no peptide links.

## Inclusion

- Pairs (i, j) with i ≠ j, `occ1[i, j] ≥ 1` and finite `q[i, j]`, wild type.
- Label: responding if `q[i, j] < 0.05`, otherwise not responding.
- No other filter. Pharyngeal and non-head neurons are kept if measured.

## Analysis plan

1. Build A_w (binary, directed; gap junctions both ways), A_p (binary,
   directed) and U = A_w OR A_p, all without self-links, on the atlas neuron
   order.
2. Compute S(A) = A + A·A for A_w and U.
3. Over the included pairs, compute AUROC (Mann–Whitney with mid-ranks) of
   S(A_w) and of S(U).
4. Split pairs by whether S(A_w) = 0 (no wired path of length ≤ 2).
5. Compute G_all, G_none and G_path as defined in CLAIM.md.
6. Cluster bootstrap: resample stimulated neurons j with replacement 2000
   times, take all included pairs of each drawn j, recompute the three gains;
   report 2.5 and 97.5 percentiles.
7. Null: 1000 degree-preserving rewirings of A_p by directed double-edge swaps
   (10 × number of edges accepted swap attempts per network, swaps that would
   create a self-link or a duplicate link are rejected). For each, rebuild U
   and recompute G_all and G_none on the same pairs and the same wired split.
8. Apply the decision rule in CLAIM.md.

Seed: 20261125 (numpy default_rng). Quick mode (CI smoke test) uses 50
bootstrap resamples and 20 null networks and does not set the status. If the
data cannot be downloaded or h5py is missing, quick mode runs the same
pipeline on a small synthetic network instead.

## Uncertainty

The atlas package gives responses pooled over animals, so a bootstrap over
animals (the unit docs/verification.md prefers) is not possible from these
files. The cluster bootstrap over stimulated neurons treats all pairs that
share a stimulated neuron as one unit. This is stated as a limit in the
README.

## Robustness checks (pre-registered)

Each is run with the same decision rule; the result table reports all of them.

- R1: the two other Ripoll-Sánchez range models.
- R2: negatives restricted to pairs that the atlas calls confidently
  non-responding (`q_eq < 0.05`); positives unchanged; pairs in neither class
  dropped.
- R3: direct-link predictor only, S(A) = A, with "no wired path" meaning no
  direct wired link.
- R4: AUROC stratified by number of trials, because pairs measured more often
  are more likely to reach q < 0.05 (Sale and Richardson 2026 report that
  observation count alone gives AUC 0.709). Strata: occurrence 1, 2, 3–4, ≥ 5.
  The stratified AUROC compares only responding and non-responding pairs in
  the same stratum: the sum over strata of the Mann–Whitney U divided by the
  sum over strata of n_responding × n_not_responding. Strata with no pair of
  one class contribute nothing. Primary network, standard labels, two-step
  predictor.

If any variant reaches the opposite pole from the primary (Supported versus
Refuted), the status is Contested. Other differences are reported but do not
change the status.

## Exploratory (does not set the status)

- unc-31 mutants are defective in dense-core-vesicle release, which carries
  neuropeptides. Among pairs with no wired path, compare the response rate of
  pairs with and without a peptide link in wild type and in unc-31.
- Counts of pairs, links and responders in each subset.

## Deviations from plan

1. Neuron order bug, fixed before the full run. The first quick run (after
   the plan was frozen) indexed matrices by the raw order of
   `neuron_ids.txt`. The atlas package actually orders neurons by
   `np.unique` of that list (sorted), which is also the order stored in
   `funatlas.h5`. A data sanity check showed the mismatch (for example AS10
   where AS2 was expected). check.py now uses the sorted order and stops if
   it differs from `funatlas.h5` other than the known AWCOF/AWCOFF spelling.
   The quick-run numbers from the buggy order were seen; nothing in the plan
   or decision rule was changed after seeing them, and a full run started
   with the bug was stopped before it printed any result.
