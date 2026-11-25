# History

## 2026-11-25: Hopfield bound, ±1 proof and novelty audit

The universality step (U) is now proved in writing ([PROOF-pm1.md](results/hopfield-rigorous-bound/PROOF-pm1.md)). That makes the energy barrier around every ±1 pattern hold for α < 0.12948, up from 0.08 on record. Its constants are certified exactly. An independent review is still pending. A novelty audit ([NOVELTY.md](results/hopfield-rigorous-bound/NOVELTY.md)) found Zhou's 2026-11-05 preprint (unrefereed). It claims ±1 fixed points up to 0.1382, which would supersede the folder's original fixed-point goal and Theorem 3. The work is retargeted at the exact energy-barrier capacity, predicted to be about 0.1298.

## 2026-11-25: Rigorous Hopfield capacity bound (Partial)

Complete proof that the Hopfield energy barrier holds for α < 0.12948 when
the interference patterns are Gaussian, with a certified constant. Matches
Stojnic's first-level value. The ±1 model is reduced to an unproved
universality hypothesis, so the rigorous ±1 bound is still 0.113.
Later the same day: Theorem 3, a reference-state refinement, proves a stable
fixed point for α < 0.134716 in the same Gaussian-interference model.

## 2026-11-25: Start the collection

Collection design, verification guide, result template and candidate list
drafted. First worked example: a finite-size simulation of Hopfield network
storage capacity (replication).
