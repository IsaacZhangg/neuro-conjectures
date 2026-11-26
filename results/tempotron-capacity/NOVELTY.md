# Novelty check

Checked 2026-11-26 (Exa page fetch, Firecrawl research paper index including
the papers citing Rubin, Monasson and Sompolinsky 2010, general web search; the
Exa search endpoint was down, so web search stood in for it). Coverage 2023 to
now plus the older literature those papers cite.

## Exact statement (the target)

Discrete-time tempotron with N inputs and T time bins. A pattern is a sequence
x(1), ..., x(T) of vectors in R^N; the neuron answers + if and only if
max_t w·x(t) > θ for some weight vector w and threshold θ. Patterns
x^μ(t), μ = 1..P, are i.i.d. standard Gaussian and the labels are i.i.d. fair
±1. The capacity α_c(T) is the load P/N below which, with probability → 1 as
N → ∞, some (w, θ) classifies all P patterns correctly, and above which none
does.

Open: any rigorous upper or lower bound on α_c(T) that tracks how it grows
with T. Counting as solved: matching bounds; counting as pushing the frontier:
new rigorous bounds, ideally with the right growth in T.

## What is known

Non-rigorous (statistical physics, simulation):
- Gütig and Sompolinsky, Nat. Neurosci. 9 (2006) 420: the tempotron; capacity
  measured in simulation, about 3 patterns per synapse at the parameters used.
- Rubin, Monasson and Sompolinsky, PRL 105 (2010) 218102, arXiv:1010.5496:
  extreme value theory plus a replica calculation on internal representations
  gives α_c ≈ ln ln K / (2 ln 2) = (1/2) log2 ln K for large K, where K counts
  effectively independent time bins (K = T in the discrete model), with a fitted
  offset; replica values for the discrete tempotron at K = 2, 3, 4. The authors
  note slow convergence and describe solution clusters. No rigorous statement.
- Baldassi, Braunstein and Zecchina, J. Stat. Mech. (2013) P12013,
  arXiv:1309.0393: tempotron-like model with discrete (binary) synapses and
  sparse binary binned inputs; replica symmetric answer α_c = 1, independent of
  K. Different synapse model, also non-rigorous.
- Yu, Tsodyks, Sompolinsky, Schmitz and Gütig, PNAS (2025),
  doi:10.1073/pnas.2426290122: tempotron with plastic short-term plasticity;
  capacity in simulation and theory for that extended model, not bounds for the
  discrete-time tempotron.

Rigorous, nearby but not this question:
- Cover (1965) counting and Gardner's (1988) calculation for the ordinary
  perceptron (T = 1): α_c = 2. Proved (Cover's function counting theorem).
- Sabato and Tishby, arXiv:1107.2021 (multiple-instance learning): worst-case
  sample complexity of learning with max over a bag of instances grows only
  poly-logarithmically in bag size. This is a VC-type statement about learning,
  not the typical-case capacity for random patterns, and gives no α_c(T).
- Gaussian comparison and capacity proofs for other perceptron variants:
  spherical and Ising perceptrons (Stojnic 2013 and later; Ding and Sun 2019;
  Huang arXiv:2404.18902; Shmalo arXiv:2607.11010), symmetric binary
  perceptrons (Aubin, Perkins and Zdeborová arXiv:1901.00314; Abbe and others),
  perceptron with variable selection (Xu, Ohzeki and Kabashima
  arXiv:2512.01861, replica). None treats a max over time bins.
- Levien and GrandPre, arXiv:2610.29531 (2026): gradient dynamics of max-pooling
  multiple-instance regression; learning dynamics, not storage capacity.

A counting upper bound in Cover's style, which uses that the tempotron's
dichotomies are determined by the signs of w·x^μ(t) − θ over all PT vectors,
gives α_c(T) ≤ about log2 T + log2 log2 T + O(1). We found it stated nowhere
for the tempotron, but it is a one-line consequence of Cover's counting and we
mark it as elementary (folklore type), not as a new result.

## What we could not rule out

- Full texts of a few theses on tempotron theory (Rubin's PhD thesis, Hebrew
  University) were not found online; they may contain sharper calculations, but
  the published papers built on them state no rigorous bounds.
- Citing papers of Rubin et al. 2010 that we only saw by abstract (learning
  rules, multi-spike tempotron, neuromorphic applications). None of the
  abstracts mentions a proof or a bound.

## What is open (and what this folder goes after)

1. A rigorous upper bound with the predicted growth (1/2) log2 ln T, well below
   the counting bound log2 T. New as far as we can find.
2. A rigorous lower bound above the trivial one, showing α_c(T) actually grows
   with T (or at least exceeds 2 for T ≥ 2). New as far as we can find.
3. Matching constants (the exact threshold). Open; we do not expect to settle
   it.

## Verdict

New as far as we can find. Every known value of the tempotron's capacity comes
from replica theory or simulation; we found no rigorous bound beyond the
elementary counting one, which we record as folklore rather than claim.
