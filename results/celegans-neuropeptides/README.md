# Do neuropeptides close the C. elegans structure–function gap?

Adding the neuropeptide connectome to the wired connectome should improve
prediction of which neuron pairs respond in the C. elegans signal-propagation
atlas, above all for pairs with no wired path.

Status: Partial
Kind: Empirical claim
Reproducibility: 🟡 runnable check with recorded output; agent review only

## Result

The peptide layer helps only where wiring is silent, and mostly through
which neurons are peptide-rich rather than which pairs are linked.

Pre-registered test on the wild-type atlas of Randi et al. 2023: 23,902
measured pairs, 1,151 responding (q < 0.05), 173 stimulated neurons. A
parameter-free path score, S(A) = A + A·A, on the wired connectome versus the
wired plus Ripoll-Sánchez short-range peptide connectome.

| Pairs | AUROC, wired | AUROC, wired + peptide | Gain [95% CI] | Degree-preserving null gain (mean) | Null p |
| :-- | :-: | :-: | :-: | :-: | :-: |
| All (G_all) | 0.608 | 0.562 | −0.046 [−0.073, −0.020] | −0.067 | 0.001 |
| No wired path ≤ 2 steps, 10,169 pairs, 314 responding (G_none) | 0.500 | 0.611 | +0.111 [+0.064, +0.157] | +0.081 | 0.001 |
| With a wired path (G_path) | | | −0.088 [−0.125, −0.051] | | |

- Over all pairs, adding the peptide layer makes prediction **worse**, so
  criterion (a) fails.
- Among pairs with no wired path, the peptide layer gives real predictive
  power (AUROC 0.61), so criterion (b) holds, and the gain is larger there
  than where a wired path exists, so (c) holds. Exactly one of (a) and (b)
  holds, so the pre-registered rule gives **Partial**.
- Most of the no-path gain survives shuffling the peptide links while
  keeping each neuron's number of peptide inputs and outputs: the shuffled
  layer already gives +0.081. Only about +0.03 AUROC is specific to which
  pairs Ripoll-Sánchez link. It beats all 1000 shuffles, but it is small.

Robustness (same decision rule):

| Variant | G_all | G_none | Null p (G_all, G_none) | Outcome |
| :-- | :-: | :-: | :-: | :-- |
| R1 mid-range peptide model | −0.098 | +0.031 [−0.010, +0.075] | 0.90, 0.23 | Refuted |
| R1 long-range peptide model | −0.109 | +0.015 [−0.028, +0.061] | 0.08, 0.19 | Refuted |
| R2 negatives only if confidently non-responding (q_eq < 0.05) | −0.036 | +0.137 | 0.001, 0.001 | Partial |
| R3 direct links only, S(A) = A | +0.003 [−0.015, +0.023] | +0.039 | 0.06, 0.007 | Partial |
| R4 AUROC stratified by number of trials | −0.045 | +0.117 | 0.001, 0.001 | Partial |

No variant reaches the opposite pole from the primary (Supported versus
Refuted), so the status stays Partial under the rule. The mid- and long-range
models, which allow peptides to travel further, carry no pair-specific
information beyond degree. Only the short-range model, which links neurons
whose processes share a nerve bundle, does. That pattern fits physical
proximity as much as peptide signalling, and this check cannot tell the two
apart.

Exploratory (does not set the status): among no-wired-path pairs, those with a
short-range peptide link respond more often than those without, 4.3% versus
2.7% in wild type and 3.9% versus 3.0% in unc-31 mutants, which lack
dense-core-vesicle release. The gap narrows in unc-31, but only a little.

## What was open and what this settles

Open: whether a neuropeptide connectome adds predictive information about the
measured signal propagation beyond wiring. Randi et al. checked 53 pairs one
by one with no null; Sale and Richardson (2026) reported a negative in one
sentence without numbers we could read. This result gives the first
pre-registered numbers with a peptide-specific null: no for the network as a
whole, a small but real yes for pairs that wiring cannot reach, and only for
the short-range model.

## Files

- [CLAIM.md](CLAIM.md): exact statement and decision rule.
- [SOURCES.md](SOURCES.md): primary sources, later work, searches.
- [METHODS.md](METHODS.md): data, analysis plan, robustness checks, deviations.
- [fetch_data.py](fetch_data.py): downloads the public data from PyPI and checks SHA-256.
- [check.py](check.py): `python check.py` (about 12 minutes on one core) or
  `python check.py --quick` (seconds; synthetic data if h5py is missing).
  Needs [requirements.txt](requirements.txt); numba speeds up the null.
- [check-output.json](check-output.json): recorded full run.
- [REVIEW.md](REVIEW.md): review record.

## Scope and limits

- One dataset: the wild-type Randi et al. 2023 atlas as bundled in
  wormneuroatlas 0.0.7.3, adult hermaphrodites, mostly head neurons. Says
  nothing about other preparations or species.
- One peptide network family (Ripoll-Sánchez et al. 2023), built from gene
  expression and receptor binding, not from measured release. The CeNGEN-based
  map in the atlas package was not tested.
- One predictor family, binary links, paths of one or two steps. A weighted or
  fitted model might use the peptide layer better.
- "No wired path" means none within two steps; almost every pair has a longer
  wired path.
- Uncertainty is a bootstrap over stimulated neurons, not over animals,
  because the bundled atlas pools animals.
- Shows association, not mechanism. The short-range result may reflect shared
  nerve bundles rather than peptides.
- One deviation from plan (a neuron-order bug fixed before the full run) is
  recorded in METHODS.md.
