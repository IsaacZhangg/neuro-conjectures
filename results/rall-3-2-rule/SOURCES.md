# Sources

Checked 2026-11-25.

## Primary

- W. Rall. "Branching dendritic trees and motoneuron membrane resistivity."
  Experimental Neurology 1:491–527 (1959). doi:10.1016/0014-4886(59)90046-9,
  PubMed 14435979. Only the abstract was read; the location of the d^(3/2)
  constraint in the full text is not verified.
- W. Rall. "Theory of physiological properties of dendrites." Annals of the
  New York Academy of Sciences 96:1071–1092 (1962). Not opened; citation from
  the standard literature.
- W. Rall. "Theoretical significance of dendritic trees for neuronal
  input-output relations" (1964), in Neural Theory and Modeling, ed.
  R. F. Reiss. Not opened; citation from the standard literature.
- W. Rall. "Rall model." Scholarpedia 4(4):1369 (2009).
  doi:10.4249/scholarpedia.1369. Read. Section "d3/2 Constraint on branch
  diameters" states the constraint, the equal-electrotonic-length condition,
  and that cat spinal motoneurons show "rough agreement" (citing Rall et al.
  1992).

## Later work that retests, refines or disputes it

- N. L. Desmond, W. B. Levy. "Dendritic caliber and the 3/2 power
  relationship of dentate granule cells." J. Comp. Neurol. 227(4):589–596
  (1984). doi:10.1002/cne.902270410, PubMed 6470225. Only the first part of the
  abstract was readable (30 Golgi-impregnated rat granule cells, mean
  caliber about 1 µm, uniform over the proximal two thirds of the molecular
  layer). The claim in docs/candidates.md that the rule holds distally but
  not proximally is **not verified**.
- A. U. Larkman, G. Major, K. J. Stratford, J. J. Jack. "Dendritic morphology
  of pyramidal neurones of the visual cortex of the rat. IV: Electrical
  geometry." J. Comp. Neurol. 323(2):137–152 (1992). doi:10.1002/cne.903230202.
  Abstract: "most dendrites did not closely obey the 3/2 branch power
  relationship".
- M. Liao, X. Liang, J. Howard. "The narrowing of dendrite branches across
  nodes follows a well-defined scaling law." PNAS 118(27):e2022395118 (2021).
  doi:10.1073/pnas.2022395118, PMC8271565. Read (results and discussion).
  Drosophila class IV neurons, superresolution diameters: per-branch-point
  exponents average 7.2 (SD 8.8), supporting no single exponent; a modified
  law d_m^p + d_0^2 = d_1^p + d_2^p fits with p ≈ 2.
- P. Desai-Chowdhry, A. B. Brummer, V. M. Savage. "How axon and dendrite
  branching are guided by time, energy, and spatial constraints." Sci. Rep.
  12:20810 (2022). doi:10.1038/s41598-022-24813-2, PMC9718790. Read (methods
  and results). NeuroMorpho.Org Golgi, Purkinje and motor neurons; mean
  dendritic daughter/parent radius ratio 0.67 ± 0.004 using section-mean
  radii, which corresponds to η ≈ 1.7 for symmetric branching (our
  conversion); Rall's value corresponds to 0.63.
- J. D. Reed, K. T. Blackwell. "Prediction of neural diameter from morphology
  to enable accurate simulation." Front. Neuroinform. 15:666695 (2021).
  doi:10.3389/fninf.2021.666695. Abstract and one passage read: many
  NeuroMorpho reconstructions lack accurate diameters; reviews reports that
  the 3/2 relation "does not always hold".
- A. D. Bird, H. Cuntz. "Optimal current transfer in dendrites." PLoS Comput.
  Biol. 12(5):e1004897 (2016). One passage read: "Rall's relationship is rarely
  satisfied in real dendrites".

What this result adds: a pre-registered test across many NeuroMorpho archives
at once, split by lab, species and cell type, with a decision rule fixed in
advance and filters for default diameters. We found no earlier study doing
exactly this, but the search below was not exhaustive.

## Datasets

- NeuroMorpho.Org (RRID:SCR_002145), REST API and CNG SWC files. License:
  CC BY 4.0 (terms of use page, read 2026-11-25). Required citations: the
  original paper for each reconstruction, NeuroMorpho.Org (RRID:SCR_002145),
  and Tecuatl, Ljungquist and Ascoli (2024), FASEB Bioadv. 6(7):207–221,
  doi:10.1096/fba.2024-00048. Download date and neuron list are recorded at
  run time.

## Searches run

- Firecrawl research index (PubMed, PMC, arXiv): "Rall 3/2 power rule
  dendritic branch points diameter test reconstructions"; "Liao Howard
  narrowing of dendrite branches across nodes scaling law exponent da Vinci
  Rall"; "Rall 1959 branching dendritic trees and motoneuron membrane
  resistivity".
- Pages read: Scholarpedia "Rall model"; PubMed 6470225; PMC8271565;
  PMC9718790 (passages); NeuroMorpho.org API reference and terms of use.
