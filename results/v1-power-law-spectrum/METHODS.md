# Methods

Committed before any code read neural responses. Before this commit, only
these were looked at: dataset metadata (experiment lists, Cre lines, imaging
depths, unit quality metrics and brain-area labels), the HDF5 node names,
shapes and dtypes of one file of each kind, and the AllenSDK source that
defines `mean_sweep_response_ns`. No response values or spike times were read.
The simulations in "Why a forward model" use synthetic data only.

## Data

Both datasets are from the Allen Brain Observatory, read from the public AWS
bucket `s3://allen-brain-observatory` (Registry of Open Data on AWS), under the
Allen Institute terms of use. Both show the same 118 natural images (plus a
blank), each for 250 ms, 50 times in random order, in back-to-back
presentations. Mice, recording method and labs' pipelines differ between them.

### D1: Visual Coding two-photon (de Vries et al., Nat. Neurosci. 2020)

- Files: `visual-coding-2p/ophys_experiment_analysis/<id>_three_session_B_analysis.h5`
  (node `analysis/mean_sweep_response_ns` and `analysis/stim_table_ns`), plus
  `visual-coding-2p/ophys_experiments.json` for metadata.
- Response per trial: Allen's `mean_sweep_response_ns`: for each presentation,
  the percent change in fluorescence relative to the mean of the preceding
  4 × 7 frames (about 0.93 s), averaged over the 14 frames (about 0.47 s) from
  stimulus onset (AllenSDK `StimulusAnalysis.get_sweep_response`, with
  `sweeplength` = presentation length, `interlength` = 4 × sweeplength,
  `extralength` = sweeplength). The running-speed column `dx` is dropped.
- Inclusion: `three_session_B` experiments whose experiment container is not
  marked failed; excitatory Cre lines only (Cux2, Emx1, Slc17a7, Rorb, Scnn1a,
  Nr5a1, Rbp4, Fezf2, Tlx3); all imaging depths. All cells in the file are
  kept, with no responsiveness filter (Stringer et al. kept all neurons).
- Unit of resampling: the mouse (`donor_id`).

### D2: Visual Coding Neuropixels (Siegle et al., Nature 2021)

- Files: `visual-coding-neuropixels/ecephys-cache/session_<id>/session_<id>.nwb`
  (nodes `units/*` and `intervals/natural_scenes_presentations/*`), plus the
  cache's `sessions.csv`, `units.csv`, `channels.csv`, `probes.csv`.
- Sessions: the 32 `brain_observatory_1.1` sessions (the functional
  connectivity sessions do not show natural scenes). One session per mouse.
- Units: AllenSDK default quality filter (`quality == good`,
  `isi_violations <= 0.5`, `amplitude_cutoff <= 0.1`, `presence_ratio >= 0.9`);
  area from the peak channel's `ecephys_structure_acronym`.
- Response per trial: spike count in [0, 250) ms from stimulus onset.
- Unit of resampling: the mouse (session).

Areas: VISp (primary), VISl, VISal, VISrl, VISpm, VISam.

Blank presentations (frame = -1) are dropped. An image with fewer than 2
trials in a session is dropped from that session's contribution; a session
missing any of the 118 images is excluded (recorded in the output).

Raw data are not committed. `check.py` downloads only the nodes it needs, by
HTTP range requests, and caches the extracted arrays outside the repository
(`$NC_CACHE`, default `~/.cache/neuro-conjectures/v1-power-law-spectrum`). The
S3 ETag, size and last-modified date of every file used are written to
`check-output.json`.

## Analysis plan

1. **Normalise.** For each neuron, divide its single-trial responses by their
   standard deviation over all 118 × 50 trials. Drop neurons with zero
   variance.
2. **Split.** For each session and each image, split its trials at random into
   two halves (sizes ⌊r/2⌋ and the rest, dropping one if r is odd so halves
   are equal), using the same split for all neurons in that session. Average
   each half: X1, X2 are 118 × N matrices for the pooled pseudo-population.
   Trial noise in neurons from different mice is independent, so pooling does
   not bias the cross-validated covariance.
3. **cvPCA.** Centre X1 and X2 over images. SVD X1 = U S Vᵀ; λ_k =
   Σ_images (X1 v_k)(X2 v_k) / 117. Repeat with the halves swapped and
   average. Average the spectrum over 10 random splits (point estimate) or
   2 splits (each bootstrap resample).
4. **Forward model (single power law).** A population spectrum
   λ_k = k^(-α), k = 1..K. Its expected 118-image cvPCA spectrum, in the large
   population limit, is the eigenvalue spectrum of Z diag(λ) Zᵀ / 117, where Z
   is a 118 × K matrix of independent standard normal stimulus features,
   centred over images. The mean log eigenvalue at each rank is computed over
   40 draws of Z (seed fixed). Ranks k ≤ 3000 are simulated exactly; the
   variance in ranks above 3000 is added as its expected value, an isotropic
   term (Σ_{k>3000} λ_k) times the centring projector (its fluctuations are
   small because it spreads over many dimensions).
   Grid: α from 0.50 to 2.00 in steps of 0.02; K ∈ {300, 1000, 3000, 10000,
   30000}.
5. **Fit.** Compare log observed λ_k with the model's mean log spectrum over
   ranks 1 to 100, with a free additive constant (overall scale), weighted by
   1 / (bootstrap variance of log λ_k + simulated variance of log λ_k); the
   bootstrap variance is computed once from the 500 resampled spectra and the
   same weights are reused for every resample. If
   some λ_k ≤ 0 among ranks 1 to 100, the fit stops at the rank before the
   first one (recorded). α̂ and K̂ are the grid point with the smallest
   weighted squared error. Each bootstrap resample is refitted the same way,
   including the choice of K.
6. **Broken power law.** Population spectrum k^(-α1) for k ≤ k_b and
   k_b^(α2-α1) k^(-α2) above, to K. Grid: α1 from 0.3 to 1.3 step 0.1, α2 from
   0.8 to 2.0 step 0.1, k_b ∈ {5, 10, 20}, K as above. Same forward model,
   same fit.
7. **Model comparison.** 100 random splits of an area's mice into two halves.
   Fit single and broken models to half A's spectrum (unweighted), compute the
   mean squared error of log λ_k over ranks 1 to 100 on half B's spectrum
   after refitting only the additive constant. The broken model is "clearly
   better" if it has the lower error in at least 95 of 100 splits and the
   median ratio of errors (broken / single) is at most 0.9.
8. **Verdicts.** Apply CLAIM.md's decision rule.

## Why a forward model

Simulations with synthetic data only (no neural data), 118 Gaussian stimuli:
a population spectrum with α = 1.04 gives a direct log-log slope of about 0.9
over ranks 1 to 10, and 0.68 to 0.98 over ranks 11 to 40 depending on how far
the power law extends (K from 10,000 down to 300). A direct slope, as in the
original paper, would therefore understate α by up to 0.4 with 118 images and
depend on an unmeasured tail. The forward model makes the finite-image effect
explicit and fits the tail length K instead of assuming it.

## Uncertainty

95% percentile intervals over 500 bootstrap resamples of mice within each area
and dataset (mice resampled with replacement; all of a resampled mouse's
neurons come along; a mouse drawn twice contributes its neurons twice). The
forward fit, including K, is redone in every resample.

## Robustness checks (planned; each is rerun for V1 in both datasets)

R1. **Fit range**: ranks 1 to 60 instead of 1 to 100.
R2. **Fewer neurons**: a random half of each mouse's neurons. The forward model
    assumes a large population; if α̂ moves by more than 0.10 the output flags
    that the large-population assumption is not met.
R3. **Response definition**: D1 with all Cre lines (including Sst, Vip,
    Pvalb); D2 with the count window [50, 300) ms.
R4. **No normalisation**: raw responses (ΔF/F percent or spike counts) without
    step 1.
R5. **Direct slope** (for comparison with the original method only, not a
    verdict): ordinary log-log slope over ranks 11 to 40 and 1 to 10.

A robustness result that moves the V1 interval outside T, or flips a dataset
between consistent and inconsistent, makes the status Contested (CLAIM.md).

## Assumptions and known differences from the original

- Allen images are 118 fixed photographs, not random draws; the forward model
  treats them as draws from an ensemble with Gaussian features.
- The forward model ignores finite neuron counts (checked by R2).
- Stringer et al. projected out 32 dimensions of spontaneous activity from
  simultaneously recorded neurons; with a pooled pseudo-population this is not
  possible, and cross-mouse trial noise is independent anyway. Shared
  within-session noise is independent across the two random trial halves in
  expectation, so it adds variance, not bias.
- Stringer et al. used deconvolved GCaMP6s activity from ~10,000 simultaneously
  recorded L2/3-4 neurons; D1 uses Allen's ΔF/F measure from GCaMP6f in
  Cre-defined populations, and D2 uses spike counts.
- 250 ms back-to-back presentations, so responses to one image overlap the
  next; image order is random, so this adds noise, not bias.

## Deviations from plan

None yet.
