# Sources

Checked 2026-11-25. arxiv.org returned 403 from the cloud environment used
for this work, so arXiv papers were read through a paper-search index
(abstracts and full-text passages) rather than from the PDFs. Journal
volume and page details marked "(not re-checked)" are from memory.

## Primary

- Feng, Shcherbina, Tirozzi. On the critical capacity of the Hopfield model.
  Comm. Math. Phys. 216 (2001), arXiv:math-ph/0012038. Statement location:
  abstract, "the probability of having a stable fixed point is one when
  α ≤ 0.113". This is the best rigorous ±1 bound we found. Its proof was
  not independently checked here.
- Stojnic. Capacity of the Hebbian-Hopfield network associative memory
  (2024), arXiv:2403.01907. Statement location: abstract. First level of
  lifting: α_c(NLT,1) = erf(x)²/(2x²) − 1 + erf(x)² ≈ 0.129490, with x fixed
  by 1 − erf(x)² − 2erf(x)e^{−x²}/(√π x) + 2e^{−2x²}/π = 0, and
  α_c(AGS,1) ≈ 0.137906. Second level (numerical): α_c(NLT,2) ≈ 0.12979,
  α_c(AGS,2) ≈ 0.138186. The analysis uses Gaussian patterns; the switch to
  binary patterns is asserted "via concentrations and the Lindeberg
  principle" without a proof in the paper.
- Newman. Memory capacity in neural network models: rigorous lower bounds.
  Neural Networks 1 (1988) (not re-checked). Energy-barrier method, α ≥ 0.056.
- Loukianova. Lower bounds on the restitution error in the Hopfield model.
  Probab. Theory Related Fields (1997) (not re-checked). α ≥ 0.071.
- Talagrand. Rigorous results for the Hopfield model with many patterns.
  Probab. Theory Related Fields 110 (1998) (not re-checked). α ≥ 0.08.
- Amit, Gutfreund, Sompolinsky. Storing infinite numbers of patterns in a
  spin-glass model of neural networks. Phys. Rev. Lett. 55 (1985). Replica
  prediction α_c ≈ 0.138 (non-rigorous).

## Later work that retests, refines or disputes it

- Steffan, Kühn (1994), Z. Phys. B (not re-checked): replica-symmetric
  corrections to the AGS value.
- Barra, Durante (2026), arXiv:2610.41938, unrefereed, proofs written with
  AI assistance: one-sided pattern universality; Corollary 9.3 compares
  Rademacher and Gaussian ground-state gains. Used in PROOF.md §6.
- Ali, Boursier, Chen, Ho, Issa, Mourrat (2026), arXiv:2611.03306,
  unrefereed: Bernoulli versus Gaussian patterns in the Hopfield free
  energy.
- Camilli, Contucci, Mingione (2026), arXiv:2610.32618, unrefereed: recent
  rigorous Hopfield work, checked for a better ±1 bound; none found.
- Genovese (2011), arXiv:1102.2535: universality of Hopfield-type free
  energies at positive temperature; does not give a zero-temperature
  capacity bound.

## Datasets

None. The result is a theorem with a certified numerical constant.

## Searches run

- Paper index: "Hopfield storage capacity rigorous lower bound", "Hopfield
  0.113 Feng Shcherbina Tirozzi", "energy barrier Newman Hopfield capacity",
  "Stojnic Hopfield fully lifted random duality", "Hopfield universality
  Bernoulli Gaussian patterns", "Hopfield model ground state Rademacher
  Gaussian comparison". Related-paper expansion from FST 2001 and
  Stojnic 2024.
- Result of the search: no rigorous ±1 bound above 0.113 was found. The
  0.1295 and 0.1379 values exist only with Gaussian patterns or with an
  unproved Gaussian switch.
