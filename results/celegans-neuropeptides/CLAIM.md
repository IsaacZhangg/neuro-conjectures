# Claim

## Statement tested

In the wild-type C. elegans signal-propagation atlas, adding the
neuropeptide connectome to the wired connectome (chemical synapses plus gap
junctions) significantly improves prediction of which stimulated → observed
neuron pairs show a response, and the improvement is largest for pairs that
have no wired path between them.

This precise version is our reading. The candidate list stated the idea in
words; no paper we found states it with a quantity (see SOURCES.md). The
closest prior statement is negative: Sale and Richardson (2026) say in one
sentence that the available neuropeptide networks "did not explain the
measured responses". This result is a pre-registered, quantitative test of
that point.

"No wired path" means no wired path of length one or two. In a nearly
strongly connected graph almost every pair has some path, so "no path at all"
would leave almost no pairs; two steps matches the predictor.

## Where it comes from

- Response data: Randi, Sharma, Dvali and Leifer (2023), the signal-propagation
  atlas. The paper reports that anatomy predicts the measured responses only
  partly and discusses extrasynaptic (peptidergic) signalling as one reason.
- Peptide network: Ripoll-Sánchez et al. (2023), the neuropeptidergic
  connectome.

Exact locations are in SOURCES.md.

## Operational definitions

- Pair: an ordered pair (j stimulated, i observed), i ≠ j, that the atlas
  measured at least once (occurrence ≥ 1 and a finite q value).
- Responding pair: q < 0.05 in the wild-type atlas (the paper's own
  threshold for a functional connection).
- Wired network A_w: a directed link j → i if the anatomical connectome
  bundled with the atlas has at least one chemical synapse j → i or at least
  one gap junction between i and j (gap junctions in both directions).
- Peptide network A_p: a directed link j → i if the Ripoll-Sánchez network
  lists at least one neuropeptide–receptor pathway from j to i.
- Combined network U: A_w OR A_p.
- Predictor (no fitted parameters): S(A) = A + A·A, the number of directed
  paths of length one or two from j to i. The predictor is scored by AUROC
  over pairs (ties count one half).
- No wired path: S(A_w) = 0 for the pair (no wired path of length one or two).
  On these pairs the wired predictor is constant, so its AUROC is exactly 0.5.

## Statistics

- G_all: AUROC of S(U) minus AUROC of S(A_w), over all pairs.
- G_none: AUROC of S(U) minus 0.5, over pairs with no wired path.
- G_path: AUROC of S(U) minus AUROC of S(A_w), over pairs with a wired path.
- Intervals: 95% percentile intervals from a cluster bootstrap over stimulated
  neurons (2000 resamples).
- Null: 1000 degree-preserving rewirings of A_p (each neuron keeps its in- and
  out-degree; no self-links or duplicate links). One-sided null p value
  p = (1 + number of null gains ≥ observed) / (1 + 1000).

## Decision rule (committed before any code reads the data)

Primary analysis: the Ripoll-Sánchez short-range model (links only between
neurons whose processes share a nerve bundle). It is the sparsest and most
anatomically constrained of the three, so its degree-preserving null is the
least saturated. We could not confirm which range model Ripoll-Sánchez et al.
treat as their main network; the mid- and long-range models are run as
robustness checks R1.

- (a) holds if the 95% interval of G_all is above 0 and its null p < 0.05.
- (b) holds if the 95% interval of G_none is above 0 and its null p < 0.05.
- (c) holds if G_none > G_path (point estimates).

- Supported if (a), (b) and (c) all hold.
- Partial if (a) and (b) hold but (c) fails, or if exactly one of (a) and (b)
  holds.
- Refuted if neither (a) nor (b) holds.
- Contested if any pre-registered robustness variant in METHODS.md reaches the
  opposite pole (Supported where the primary is Refuted, or the reverse).
