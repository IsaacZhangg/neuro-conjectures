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

- A. Wendemuth, "Storage capacity bounds in multilayer neural networks",
  Int. J. Neural Systems 5 (1994) 217–228. Its abstract says it gives "general
  analytic expressions for the lower and upper storage capacity bounds" for
  networks with a fixed output function, including committee machines, and that
  simulations sit "just slightly above the lower limit". The full text is
  paywalled and we could not read it, so we cannot tell whether its lower limit
  is rigorous or what value it gives for K = 3. It is from a physics group and
  cited once (West and Saad 1998, a statistical mechanics paper), which suggests
  a heuristic bound, but this is unverified. Checking it is the first open item.

## What is open

- The exact value of α_c(3) (gap now 1.55 to 3.43).
- Any proved lower bound for ±1 inputs.
- Proved lower bounds closer to the physics prediction (about 3).

## Verdict

**New, as far as we can find**, with one unread 1994 paper as the caveat. The
result is a first proved lower bound, not a new value: physics already predicts
a capacity near 3 and simulations store about 2, so 1.55 confirms rigorously
part of what was expected rather than discovering something unexpected. The
natural next target is the open part: raising the proved lower bound toward 2
or beyond, and the ±1 case.
