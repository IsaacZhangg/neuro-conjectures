# Sources

Checked 2026-11-25. "Read" means the full text or the named passages were
opened; "snippet" means only an abstract or search excerpt was seen.

## Primary

- Cuntz H, Mathy A, Häusser M. A scaling law derived from optimal dendritic
  wiring. PNAS 109(27):11014–11018, 2012. doi:10.1073/pnas.1200430109. Read
  (PMC3390826). Statement location: Eq. 1 and its derivation in Results; the
  branch-point exponent 0.72 ± 0.10 (74 NeuroMorpho datasets with ≥ 10 cells,
  0.70 pooled) in the NeuroMorpho analysis. Not verified: the equation image
  itself (the image server was refused by this environment's network proxy)
  and supplementary Fig. S7 (2D case).

## Later work that retests, refines or disputes it

- Ferreira Castro A, Cardona A. Invariant synaptic density across species.
  bioRxiv 2024, doi:10.1101/2024.07.18.604056; eLife Reviewed Preprint
  RP108642, doi:10.7554/eLife.108642.1. Read. **Prior EM test**: compares the
  length predicted by Eq. 1 (n = dendritic synapses, fixed exponents) with
  measured length for 2,486 EM neurons (fly larva, adult fly hemibrain plus
  ellipsoid body, zebrafish olfactory bulb, mouse MICrONS); correlations
  0.64 to 0.95, with the law overestimating small and underestimating large
  arbours. No free-exponent or per-cell-type fits; FlyWire excluded as too
  incomplete. This result differs in testing the branch-point form, fitting
  the exponent, covering FlyWire whole brain, and checking sensitivity to
  pruning and to n.
- Anton-Sanchez L et al. PLoS Comput Biol 2018, doi:10.1371/journal.pcbi.1006593.
  Read (passages). Uses EM fly data for a regularity index; its scaling
  statement comes from models, not from fitting EM data.
- Otopalik AG et al. eLife 2017, doi:10.7554/eLife.22352. Read (passages).
  Crab stomatogastric neurons are highly tortuous and, by the authors'
  reading, do not follow wiring optimisation (light microscopy).
- Cuntz H et al. PLoS One 2013, doi:10.1371/journal.pone.0071540. Snippet.
  Fly HS cells follow the planar (1/2) law.
- Ferreira Castro A et al. eLife 2020, doi:10.7554/eLife.60920. Read. Fly
  sensory dendrites follow the 2D law through development.
- Lu J et al. PLoS Biol 2009, doi:10.1371/journal.pbio.1000032. Snippet.
  Motor axons consistent with the 2D law.
- Bird AD, Cuntz H. Cell Rep 2019, doi:10.1016/j.celrep.2019.04.097. Read
  (passages). Builds on the law; no test.
- Fritz et al. bioRxiv 2025, doi:10.1101/2025.07.26.666968. Snippet only; not
  verified whether it tests the law.

## Datasets

- FlyWire whole-brain connectome, materialization 783: skeletons in the public
  bucket `gs://flywire_v141_m783/skeletons_mip_1` (read via
  storage.googleapis.com). Dorkenwald et al., Nature 2024,
  doi:10.1038/s41586-024-07558-y. License: not verified in this session (the
  FlyWire site was unreachable from here); FlyWire states its data are openly
  available for reuse with citation.
- FlyWire annotations: github.com/flyconnectome/flywire_annotations, commit
  a83b2776d60d5764cef36b927f5f9679c16c47a2. Schlegel et al., Nature 2024,
  doi:10.1038/s41586-024-07686-5; cite also Matsliah et al. 2024 and Berg et
  al. 2026 per that repository. License: no license file in the repository;
  not verified.
- MICrONS minnie65 (MICrONS Consortium, Nature 2025): skeletons v661 and cell
  types in `s3://bossdb-open-data/iarpa_microns/minnie/minnie65/`. License:
  not verified in this session (microns-explorer.org was unreachable); the
  MICrONS release is described by its authors as CC BY 4.0.

## Searches run

PubMed/PMC full text and citing papers of doi:10.1073/pnas.1200430109 via a
paper index; web searches for "2/3 power law dendrite length branch points
connectome", "wiring optimization EM reconstruction scaling", "Cuntz scaling
law hemibrain / FlyWire / MICrONS / H01", and Cuntz-lab papers 2013–2025. The
paper index listed about 60 citing papers against about 125 reported by
Google Scholar, so the search is not exhaustive.
