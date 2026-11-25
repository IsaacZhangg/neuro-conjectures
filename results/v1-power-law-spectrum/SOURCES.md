# Sources

Checked 2026-11-25.

## Primary

- Stringer C, Pachitariu M, Steinmetz N, Carandini M, Harris KD.
  High-dimensional geometry of population responses in visual cortex.
  Nature 571, 361-365 (2019). doi:10.1038/s41586-019-1346-5. Open text:
  PMC6642054. Statement location: Fig. 2c ("1/n^α where α = 1.04"),
  Fig. 2e (1.04 ± 0.017 SE, n = 7 recordings), Methods "Estimation of power-law
  exponent" (fit over ranks 11 to 500). Low-dimensional stimuli: Fig. 3e-g
  (α = 1.49 for 8-d, 1.65 for 4-d projections, 3.43 for gratings; bound
  α > 1 + 2/d). Read in the bioRxiv version, doi:10.1101/374090.

## Later work that retests, refines or disputes it

- Pospisil DA, Pillow JW. Revisiting the high-dimensional geometry of
  population responses in the visual cortex. PNAS (2025).
  doi:10.1073/pnas.2506535122 (preprint doi:10.1101/2024.02.16.580726).
  Reanalyses the Stringer data; argues cvPCA is biased and a broken power law
  (slope about 0.5 then about 1.2) fits better; average tail α 1.20 with their
  moment estimator versus 1.01 with cvPCA. Read in the preprint; the
  breakpoint value and comparison criterion were not read verbatim.
- Davidovich IA, Roudi Y. Bayesian interpolation for power laws in neural data
  analysis. arXiv:2204.08525 (2022). Fitted exponents on the Stringer data
  depend on the noise model, fit range and recording. Abstract only.
- Kong NCL, Margalit E, Gardner JL, Norcia AM. Increasing neural network
  robustness improves match to macaque V1 eigenspectrum... PLoS Comput Biol
  (2022). doi:10.1371/journal.pcbi.1009739. Macaque V1 cvPCA exponent near 1.
  Abstract only.
- Gauthaman RM, Ménard B, Bonner MF. Universal scale-free representations in
  human visual cortex. PLoS Comput Biol (2025). doi:10.1371/journal.pcbi.1013714.
  Human fMRI (Natural Scenes Dataset), power-law variance decay across visual
  regions. Abstract only.
- Ghosh et al. Why all roads don't lead to Rome: representation geometry varies
  across the human visual cortical hierarchy. arXiv:2509.13459 (2025). Human
  NSD fMRI; α close to 1 in areas with a power-law fit, but some anterior areas
  do not fit a power law.
- Pachitariu et al. 2026, Nature, doi:10.1038/s41586-026-10528-1, reported by
  the conjecture-finder channel to find CA1 responses closer to uncorrelated.
  Not opened in this session.

## Datasets

- D1: Allen Brain Observatory Visual Coding two-photon. de Vries SEJ et al.,
  Nat. Neurosci. 23, 138-151 (2020). doi:10.1038/s41593-019-0550-9.
  `s3://allen-brain-observatory/visual-coding-2p/` (files last modified
  2019-10; exact ETags in check-output.json). Allen Institute terms of use.
- D2: Allen Brain Observatory Visual Coding Neuropixels. Siegle JH et al.,
  Nature 592, 86-92 (2021). doi:10.1038/s41586-020-03171-x.
  `s3://allen-brain-observatory/visual-coding-neuropixels/ecephys-cache/`.
  Allen Institute terms of use.
- Not reachable from this environment (proxy refused): the original Stringer
  data (Janelia figshare, doi:10.25378/janelia.6845348), Zenodo, the DANDI
  API, OSF, GIN, CRCNS. The Natural Scenes Dataset needs a data use agreement
  and was not used.

## Searches run

- Firecrawl research search: "eigenspectrum power law visual cortex population
  responses natural images alpha exponent broken power law cross-validated PCA"
  (from 2019-06).
- Web search: "Pospisil Pillow broken power law V1 eigenspectrum Stringer data";
  "cvPCA eigenspectrum power law Allen Brain Observatory natural scenes
  two-photon alpha exponent"; "DANDI dandiset Stringer Pachitariu natural
  images 2800".
- DANDI metadata (all 1,294 dandisets, via the public S3 bucket), searched for
  natural-image, Stringer/Pachitariu, ImageNet and mesoscope keywords. No
  reachable dandiset with many natural images and repeats in mouse cortex was
  found besides Allen's.
- No published cvPCA exponent for the Allen natural-scene sessions was found.
