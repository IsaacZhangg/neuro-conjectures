# Candidate targets

Screened targets that nobody has started. Each needs a full screening issue
(see DESIGN.md) before work begins. Sources were checked on 2026-11-25 from the
pages linked; statement locations still need confirming against the full text.

## 1. Rall's 3/2 power rule at dendritic branch points

In progress: plan frozen in [results/rall-3-2-rule/](../results/rall-3-2-rule/)
(branch `feat/rall-3-2-rule`).

- **Claim.** At a branch point, the parent diameter to the 3/2 power equals the
  sum of the daughter diameters to the 3/2 power. This is the condition that
  lets a branched tree be collapsed into an equivalent cylinder.
- **Kind.** Empirical claim (the rule is a modelling condition; the question is
  whether real dendrites satisfy it).
- **Sources.** Rall (1959, 1962, 1964). Summary in the Scholarpedia article
  "Rall model". Partial tests: Desmond and Levy (1984), dentate granule cells,
  rule holds distally but not proximally (PubMed 6470225); Mainen et al.,
  CA1 pyramidal reconstructions.
- **Test.** Branch-point exponent fits across NeuroMorpho.org reconstructions,
  split by cell type, species and lab, restricted to reconstructions with
  measured (not default) diameters. Estimate the exponent per cell and test
  whether 1.5 lies in the interval.
- **Why it matters.** Equivalent-cylinder reductions are used widely in compact
  neuron models. A systematic deviation tells modellers where the reduction is
  unsafe.
- **Data.** NeuroMorpho.org (open, per-archive licenses).

## 2. Power-law eigenspectrum of visual cortex responses (α ≈ 1)

- **Claim.** The variance of the n-th principal component of mouse V1 responses
  to natural images decays as 1/n^α with α just above 1 (reported 1.04 ± 0.02 SE,
  7 recordings).
- **Kind.** Empirical claim, currently Contested in the literature.
- **Sources.** Stringer, Pachitariu, Steinmetz, Carandini and Harris, Nature
  2019, doi:10.1038/s41586-019-1346-5 (open text PMC6642054). Pospisil et al.,
  PNAS 2025, doi:10.1073/pnas.2506535122, argue a broken power law
  fits better.
- **Test.** Apply the cross-validated PCA method to an independent dataset,
  for example Allen Brain Observatory Visual Coding natural scenes sessions,
  and compare single and broken power-law fits with a pre-registered criterion.
- **Why it matters.** The exponent bounds how smooth the population code can be;
  the paper argues that α > 1 + 2/d is needed for a smooth (differentiable) code of d-dimensional stimuli.
- **Data.** Allen Brain Observatory (open), original data from the authors.

## 3. Grid cell module scale ratio

- **Claim.** Adjacent grid modules form a geometric progression with ratio
  √e ≈ 1.65 under winner-take-all decoding, or about 1.44 under probabilistic
  decoding, with modest variability within and between animals.
- **Kind.** Empirical test of a theoretical prediction.
- **Sources.** Wei, Prentice and Balasubramanian, eLife 2015, 4:e08362
  (PMC4616244). Comparison data: Barry et al. 2007 (mean 1.64, n = 6);
  Stensola et al. 2012 (mean 1.42 ± 0.17, n = 24).
- **Test.** Recompute module periods and ratios from openly shared grid cell
  recordings with several modules per animal, and test constancy of the ratio
  within animals.
- **Why it matters.** Distinguishes between proposed decoding schemes for the
  grid code.
- **Data.** To be found during screening (candidates: CRCNS hc and mec sets,
  DANDI grid cell dandisets). This is the main risk for this target.

## 4. Universal scaling of cortical folding

- **Claim.** Across mammals, total cortical area A_t, exposed area A_e and
  thickness T obey A_t · T^(1/2) = k · A_e^(5/4), with a single universal
  exponent of about 1.25.
- **Kind.** Empirical claim, disputed in published comments.
- **Sources.** Mota and Herculano-Houzel, Science 2015, 349:74 (PubMed
  26138976); comment by Lewitus et al. 2016 and the authors' response,
  doi:10.1126/science.aad2346 (which reports an exponent of 1.305 on extended
  data).
- **Test.** Refit within humans using open structural MRI (for example the
  Human Connectome Project or OpenNeuro datasets with FreeSurfer outputs), and
  test whether the within-species exponent matches the across-species value.
- **Why it matters.** Tests whether folding follows a simple physical rule.
- **Data.** HCP requires accepting data use terms; OpenNeuro datasets are CC0.

## 5. Neuronal avalanche size exponent of 3/2

- **Claim.** Spontaneous cortical activity organises into avalanches whose sizes
  follow a power law with exponent about −3/2, as expected at a critical
  branching point.
- **Kind.** Empirical claim, Contested.
- **Sources.** Beggs and Plenz, J. Neurosci. 2003, 23:11167. Later debate on
  subsampling, binning and alternative explanations (to be compiled during
  screening).
- **Test.** Pre-registered avalanche analysis on open spiking data with several
  bin widths and a likelihood-ratio comparison against lognormal and
  exponential fits.
- **Data.** DANDI or CRCNS spiking datasets.

Note: sources for target 5 have not been opened in this session; confirm them
before screening.
