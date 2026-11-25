# Novelty check

Checked 2026-11-25 (Exa web search, Firecrawl research paper index with citing
papers of arXiv:2312.08244 and arXiv:2312.08257, general web search). Coverage
2023 to now plus the older literature the recent papers cite.

## Exact statement

For the K = 3 tree committee machine (three sign perceptrons on disjoint thirds
of n inputs, majority output) with i.i.d. Gaussian patterns and arbitrary labels,
all ⌊αn⌋ patterns are storable with probability → 1 for every α ≤ 1.55. So
α_c(3) ≥ 1.55.

## What is known

Rigorous, upper side:
- Mitchison and Durbin (1989): α_c(3) ≤ about 5.42.
- Stojnic, arXiv:2312.08244 (2023): ≤ 4.025, equal to the replica symmetric value.
- Stojnic, arXiv:2312.08257 (2023; Information and Inference 15(3), 2026,
  doi:10.1093/imaiai/iaag034): ≤ 3.43 for K = 3. Best proved upper bound found.
- Stojnic's later papers (arXiv:2402.05696, arXiv:2402.05719 / Neural Networks
  2026) treat other activations or wide networks. Upper bounds only.

Non-rigorous (statistical physics, simulation):
- Replica symmetric 4.02 (Barkai, Hansel and Sompolinsky 1992; Engel et al.
  1992); one-step replica symmetry breaking about 3.0 to 3.14 (same, and Malzahn
  and Engel, arXiv:cond-mat/9906201).
- Monasson and Zecchina (arXiv:cond-mat/9501082): replica symmetry argued
  exact below α ≈ 1.8 for K = 3; internal-representation entropies vanish at
  2.9 and 3.5.
- Simulations store about 2 to 2.75 patterns per weight (Engel et al. 1992:
  "slightly more than twice the number of degrees of freedom").

Lower side, rigorous: none found. Both Stojnic papers list only upper bounds.
Scaling-type memorization results for large networks (Baum 1988 and later work
on threshold and ReLU networks) use other architectures with free output
weights or thresholds and give no constant for this model.

## What we could not rule out

A. Wendemuth, "Storage capacity bounds in multilayer neural networks",
Int. J. Neural Systems 5(3) (1994) 217–228, doi:10.1142/S0129065794000232
(also a talk at EUFIT 94, and his 1994 Oxford thesis "Optimisation in neural
networks").

- Abstract: "General analytic expressions are given for the lower and upper
  storage capacity bounds of multilayer neural networks which have variable
  weights between input layer and first hidden layer, and a fixed output
  function ... The special cases of committee and parity machines ... are
  discussed ... An explanation is given as to why [simulations] have a storage
  capacity just slightly above the lower limit and how this can be improved."
- Full text: paywalled. No open preprint found (searched Exa and the web for
  the paper, the conference version and the thesis; the thesis is listed in
  EThOS/OpenGrey with no open copy).
- Citing papers: Crossref lists one, West and Saad, J. Phys. A 31 (1998) 8977
  (paywalled). West's open 1997 Edinburgh PhD thesis ("Role of biases in neural
  network models", era.ed.ac.uk handle 1842/11546), chapter 4, cites
  Wendemuth's bounds work together with Baum and Haussler 1989, Mitchison and
  Durbin 1989, Bartlett 1993, Maass 1994 and Sakurai 1995, and sums up that
  literature as "the well known lower and upper bounds of the storage capacity
  (per adjustable weight) of 1 and log₂ K". No K = 3 committee value is quoted.
  Monasson–Zecchina and Engel–Van den Broeck-era papers we read do not quote
  Wendemuth's bound.
- So the only lower bound we can attribute to this literature is 1 per
  adjustable weight, below 1.55. A fixed-internal-representation counting
  argument with Cover's theorem gives exactly 1 for K = 3, which fits that
  summary, but we could not confirm that this is what Wendemuth computed, nor
  whether his "lower limit" for K = 3 is rigorous or larger. **Not checked:**
  if the paper gives a rigorous K = 3 bound of 1.55 or more, this result is
  partly known.

## What is open

- The exact value of α_c(3) (gap now 1.55 to 3.43).
- Any proved lower bound for ±1 inputs.
- Proved lower bounds closer to the physics prediction (about 3).

## Verdict

**New, as far as we can find.** The only lower bound we can trace in the older
literature is 1 per weight; one 1994 paper (Wendemuth) could not be read and is
the remaining caveat. The
result is a first proved lower bound, not a new value: physics already predicts
a capacity near 3 and simulations store about 2, so 1.55 confirms rigorously
part of what was expected rather than discovering something unexpected. The
natural next target is the open part: raising the proved lower bound toward 2
or beyond, and the ±1 case.
