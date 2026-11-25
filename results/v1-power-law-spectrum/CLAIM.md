# Claim

## Statement tested

In mouse visual cortex, the stimulus-related variance of the n-th principal
component of a large population's responses to natural images falls off as a
power law, λ_n ∝ n^(-α), with α just above 1 (Stringer et al. report
α = 1.04 ± 0.02 SE over 7 recordings of about 10,000 V1 neurons and 2,800
images).

## Where it comes from

Stringer, Pachitariu, Steinmetz, Carandini and Harris, "High-dimensional
geometry of population responses in visual cortex", Nature 571, 361-365 (2019),
doi:10.1038/s41586-019-1346-5. Fig. 2c (eigenspectrum, "1/n^α where
α = 1.04"), Fig. 2e (1.04 ± 0.017 SE, n = 7 recordings). Exponent fitted by
linear regression of log variance on log rank over ranks 11 to 500
(Methods, "Estimation of power-law exponent"). Locations confirmed against the
bioRxiv version (doi:10.1101/374090), not the journal PDF.

## Our reading (stated because the original data are not reachable here)

The original recordings (Janelia figshare 6845348) cannot be downloaded from
this environment. The only open natural-image recordings reachable here are
the Allen Brain Observatory sessions, which show 118 natural images 50 times
each. With 118 images a spectrum has at most 117 ranks, and simulations run
before any data were read (see METHODS.md, "Why a forward model") show that a
direct log-log slope over ranks 11 to 40 of a 118-image spectrum is far below
the population exponent (about 0.6 to 0.8 when the true α is 1.04). A direct
slope therefore cannot test the claim.

So the precise statement tested here is:

> The cvPCA eigenspectrum of mouse visual cortex responses to the 118 Allen
> natural images is consistent with 118 random draws from a stimulus ensemble
> whose population spectrum is a single power law n^(-α) with α in
> [1.00, 1.10], and a broken (two-slope) power law does not fit clearly better.

This is our reading. It replaces "the spectrum is a power law over ranks
11-500" with "the finite-sample spectrum matches what such a power law would
produce", under the forward-model assumptions listed in METHODS.md.

## Operational definitions

- **Population response**: responses of all included neurons in one visual
  area, pooled across mice into one pseudo-population (all mice saw the same
  118 images).
- **cvPCA spectrum**: for each neuron, the 50 trials of each image are split
  into two halves and averaged; principal directions are found from one half
  and the variance along them is measured as the covariance with the other
  half (Stringer et al. 2019, Methods), symmetrised over the two halves.
- **α̂**: the exponent of the single power-law population spectrum whose
  simulated 118-image cvPCA spectrum best matches the observed one over ranks
  1 to 100 (forward model in METHODS.md), with the tail length K chosen jointly.
- **Interval**: 95% percentile interval over 500 bootstrap resamples of mice.
- **Broken power law clearly better**: across 100 random splits of mice into
  halves, the broken-power-law forward model fitted to one half predicts the
  other half's spectrum with lower error in at least 95% of splits, and the
  median error is at least 10% lower.

## Decision rule (committed before any code reads the data)

Primary area: VISp (V1). Two datasets: Allen Brain Observatory Visual Coding
two-photon (D1) and Visual Coding Neuropixels (D2).

Bands: claim band B = [1.00, 1.10]; tolerance band T = [0.90, 1.20].

Per dataset, V1 is:

- **consistent** if the 95% interval of α̂ lies inside T and overlaps B, and
  the broken power law is not clearly better;
- **inconsistent** if the 95% interval of α̂ lies entirely outside T;
- **inconclusive** otherwise (including: broken power law clearly better).

Overall status:

- **Partial** (the strongest status this test can give, because 118 images
  test a narrower scope than 2,800) if V1 is consistent in both datasets and
  no robustness check in METHODS.md moves either interval outside T.
- **Refuted** if V1 is inconsistent in both datasets and this survives the
  robustness checks.
- **Contested** if the two datasets disagree (one consistent, one
  inconsistent), if the broken power law is clearly better in V1 in both
  datasets, or if a robustness check flips a dataset between consistent and
  inconsistent.
- Otherwise **Partial**, labelled inconclusive, with the reason.

Secondary (reported, does not set the status): the same per-area verdict for
VISl, VISal, VISrl, VISpm and VISam in both datasets, as a test of whether
α ≈ 1 holds across mouse visual areas.
