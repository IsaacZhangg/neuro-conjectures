# History

## 2026-11-26: Hopfield bound, ±1 proof passes agent review

A separate agent reviewed the ±1 energy-barrier proof ([PROOF-pm1.md](results/hopfield-rigorous-bound/PROOF-pm1.md)) and its certificate independently and found no gap. Six cosmetic issues were fixed. The status stays Partial until a person checks the proof.

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

## 2026-11-25: Neuropeptides and C. elegans signal propagation

Pre-registered test of whether the Ripoll-Sánchez neuropeptide connectome
improves prediction of the Randi et al. 2023 signal-propagation atlas beyond
wiring. Status Partial: it lowers AUROC over all pairs, raises it to 0.61 for
pairs with no wired path within two steps, and most of that gain comes from
peptide degree. Agent-drafted and agent-reviewed only.

## 2026-11-25: Dense associative memory capacity, plan frozen

Claim and proof plan for the all-pattern capacity of dense associative memory
with Gaussian patterns (Lucibello and Mézard 2024 open question) committed
before proof work. Status Open.

## 2026-11-25: Dense associative memory capacity, proved

Proof that the all-pattern capacity of dense associative memory with Gaussian
patterns equals the Lucibello and Mézard union bound for every λ ∈ (0, ∞]
under the static retrieval criterion, with constants certified by Arb and an
independent AI-agent review. Status Proved.

## 2026-11-25: SCS network global stability (Partial)

Proved that the silent state of the Sompolinsky–Crisanti–Sommers random
network is globally stable with high probability for g < 1/√2, past the
known g < 1/2. The range 1/√2 ≤ g < 1 stays open. Not yet independently
reviewed.

## 2026-11-25: Dense associative memory capacity, novelty check

Added NOVELTY.md. Fresh literature search found no earlier rigorous proof
that the Gaussian all-pattern capacity equals the union bound; the union
bound itself is Lucibello and Mézard's and is credited to them. Verdict: new.

## 2026-11-25: C. elegans neuropeptides, novelty check

Added NOVELTY.md. The wired path-score baseline with a degree-preserving null
and the one-line negative for peptide networks are already in Sale and
Richardson (2026); the pre-registered peptide-specific test, the
no-wired-path gain and the proximity control are not, as far as we found.
Verdict: partly known.

## 2026-11-25: SCS global stability, novelty audit

Added results/scs-global-stability/NOVELTY.md. No rigorous bound above g = 1/2 for the SCS model was found in print, but the same argument and constant 1/√2 are published for random Lotka–Volterra systems (Clenet, Massol and Najim 2023). Verdict: partly known. Status unchanged (Partial).

## 2026-11-25: Tree committee machine, K = 3

New result: a proved lower bound α_c(3) ≥ 1.55 for the storage capacity of the
K = 3 tree committee machine with Gaussian patterns, the first explicit proved
lower bound we found. Status Partial: the exact capacity is still open (best
proved upper bound 3.43). Agent-drafted, review pending. NOVELTY.md
verdict: new as far as we can find; Wendemuth (1994) could not be read.
