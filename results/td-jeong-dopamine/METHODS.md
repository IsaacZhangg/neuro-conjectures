# Methods

Frozen 2026-11-25, before any code read the photometry or event data. Only the
DANDI file listing (asset names and sizes) had been seen.

## Data

DANDI 000351, version 0.260828.1905 (CC-BY-4.0), "Mesolimbic dopamine release
conveys causal associations". NWB files per session, each with an event log
(event code, time in ms, no-solenoid flag) and processed dF/F
(`processing/photometry/dff`).

Animals and the tests they enter, as in Jeong's analysis code
(github.com/namboodirilab/ANCCR, `analysis/`):

- HJ-FP-F1, F2, M2, M3, M4, M6, M7, M8: RandomRewards sessions (Tests 1, 2).
  M8 has RandomRewards only.
- HJ-FP-M2, M3, M4, F1, F2, M6, M7 (this order): Pavlovian, 8 s cue,
  background reward, extinction and Poisson-cue sessions (Tests 3 to 9). Session
  ranges are Jeong's `startday`, `endday`, `learnedday` tables in
  `fig4/cumsum_pavlovian.m`, `fig4/cumsum_bgdrw_vs_extinction.m`,
  `figS10/omission.m`, `fig5/intermediate_rsp_trialless.m`.
- Sequential conditioning cohort, controls HJ-FP-datWT-stGtACR-F1 and
  HJ-FP-WT-stGtACR-F1, F2, F3, M1, M2, M3 (Test 10), and inhibited animals
  HJ-FP-datHT-stGtACR-M2 to M7 (Test 11), as in `fig6/*sequential*.py`.

Only the event log and dF/F are read. Raw photometry, which makes the opto files
large, is not downloaded: those files are read by HTTP range requests. No raw
data is committed. `fetch_data.py` (to be added) records file names, sizes and
SHA-256 sums of what was read.

Event codes used: 5 lick onset, 7 background (uncued) solenoid, 10 cued
solenoid, 15 sound 1, 16 sound 2, 0 session end; laser codes are taken from each
file's `eventidx_table`. Reward time is the first lick within 5 s after a
solenoid (as in Jeong's code); if there is none, the solenoid time. Cued
solenoid entries flagged "no solenoid" are omissions. CS+ is the sound that is
followed by reward on the first Pavlovian session (Jeong's rule). In the
sequential task CS1 and CS2 are the odd and even code-15 events (Jeong's rule);
if the event table gives them different codes we use those.

## Model family TD-MS

Time is discretised at Δ = 0.1 s. Every event in the log is placed in the bin
that contains it. Within a session the model sees every event in order;
between sessions all traces (features, eligibility) reset to zero and weights
carry over.

Event types: CS+ onset, CS+ offset, CS− onset, CS− offset, reward (cued and
uncued sucrose are one type, since they are the same sucrose), and in the
sequential task CS1 onset/offset and CS2 onset/offset. Licks are actions and are
not features. Offset events exist only when the cue has a duration (sound off is
CS onset + CSduration from the log).

Features. For each event type e and microstimulus i = 1..m,

    x_{e,i}(t) = Σ_k  y(t − t_k) · exp(−(y(t − t_k) − μ_i)² / (2σ²)),   y(τ) = d^τ,

summed over all past occurrences t_k ≤ t of e within the session (linear
superposition, cut at 60 s), μ_i = i/m, d the memory decay per second, σ the
width (Ludvig, Sutton & Kehoe 2008). Optional features: a constant bias
feature (b ∈ {0, 1}); a complete-serial-compound representation (one indicator
per 0.1 s bin since each event, up to 10 s, superposed) in place of
microstimuli, as the second representation class.

Value and learning, with r_t = 1 in the bin of a reward:

    V_t = w · x_t
    δ_t = r_t + γ V_t − V_{t−1} − I · L_t
    e_t = γ λ e_{t−1} + x_{t−1}
    w  ← w + α g(δ_t) e_t,   g(δ) = δ if δ ≥ 0 else η δ.

γ = γ_s^Δ with γ_s the discount per second. L_t = 1 while the laser is on (only
the inhibited sequential-task animals), and I is fixed per animal by calibration
(below), not fitted. Initial weights all equal w₀ (w₀ > 0 is optimistic
initialisation).

Dopamine. Simulated dF/F = s · (δ ∗ h)(t), with h(τ) = e^{−τ/0.5 s} − e^{−τ/0.05 s}
scaled to peak 1. The scale s is per animal: the mean data reward response
(AUC 0 to 0.5 s after reward minus the pre-cue baseline) over that animal's
last three sessions before `learnedday` divided by the same model measure. For
the sequential cohort the reference is their first Pavlovian session.

Measurement noise. Every per-trial measure m (a cue AUC, a reward AUC) gets
independent Gaussian noise with SD equal to that animal's data SD of
(m − running median of m over 21 trials) in the same sessions. Each test is run
with 20 noise seeds (seeds 0 to 19); the reported t and p are the medians over
seeds. Without this, a deterministic model gives AIC and t values with no
measurement error and the tests would not mean the same thing as on data.

Inhibition calibration (Test 11). For each inhibited animal, I is set so that the
model's ratio of (AUC 0 to 0.5 s from CS2 + 1 s, minus baseline) to the reward
response on the first Pavlovian session equals that animal's data ratio (Jeong's
`measure_inhibition_sequential.py`, about −0.6 on average).

Free parameters and search box (the same box for every animal):

| Parameter | Range |
| :-- | :-- |
| α learning rate | 1e−3 to 0.3, log scale |
| γ_s discount per second | 0.5 to 0.999 |
| λ trace | 0 to 0.99 |
| m microstimuli | 5, 10, 20 |
| σ width | 0.02 to 0.2 |
| d memory decay per second | 0.3 to 0.99 |
| w₀ initial weight | 0 to 0.5 |
| b bias feature | 0 or 1 |
| offset features | on or off |
| η negative-error scale | 0.02 to 1 |
| representation | microstimulus or CSC |

Fixed and not tuned: Δ, the kernel h, the 60 s superposition cut, reward size,
the noise procedure, the scale rule, and all analysis windows.

## Tests and pass rules

Each test reuses Jeong's measure, window and statistic from the named script.
"Same sign, p < 0.05" means the model's t has the sign of the paper's t and the
two-tailed p is below 0.05 on the same n. Where the paper reports a null,
"p ≥ 0.05" is the pass. For Tests 9 to 11, where the paper gives no statistic,
the pass rule uses the data's own 95% confidence interval across animals
(t-based), computed by our code from the same files.

1. **Sucrose response rises over RandomRewards sessions** (`fig3/corr_randomrewards.m`).
   Response = AUC −0.5 to 1 s minus AUC −2 to −0.5 s around the first lick after
   each background solenoid, excluding IRI < 3 s and rewards without a lick, 12 s
   IRI sessions only. Per animal, correlation with reward number across
   sessions; one-sample t, n = 8. Pass: t > 0, p < 0.05 (paper t(7) = 4.40).
2. **Correlation with previous IRI** (same script). Per session correlation of
   response with the previous IRI, averaged per animal; one-sample t, n = 8.
   Pass: t > 0, p < 0.05 (paper t(7) = 5.95).
3. **Dopamine learning is more gradual than behaviour** (`fig4/cumsum_pavlovian.m`,
   condition 1). Model DA cue measure (AUC 0 to 2 s / 2 minus AUC −1 to 0 s)
   against the animal's real lick measure; normalised, smoothed cumsums;
   abruptness = max distance from the chord; change trial = its argmax. Paired t,
   n = 7. Pass: abruptness same sign as paper (t(6) = 9.06), p < 0.05, and change
   trial same sign as paper (t(6) = −2.93), p < 0.05.
4. **No dopamine change when the cue is lengthened to 8 s** (same script,
   condition 2). One-sample t of the model DA ΔAIC (constant minus Weibull fit,
   Jeong's bounds). Pass: p ≥ 0.05 or t ≤ 0 (paper t(6) = −0.86, p = 0.42).
5. **No dip at the old delay, dip at omission** (`figS10/omission.m`). Cumulated
   model response AUC 3 to 4 s after CS+ in the first 8 s session versus AUC 9 to
   10 s in the first extinction session, Jeong's baselines. Pass: paired t has
   the sign of the data's paired t with p < 0.05, and the 3 to 4 s response is not
   significantly negative (one-sample t ≥ 0 or p ≥ 0.05).
6. **Cue response stays positive through abrupt extinction of licking**
   (`fig4/cumsum_pavlovian.m`, condition 3). Pass: abruptness paired t same sign
   as paper (t(6) = 5.67), p < 0.05, and the mean model cue response over the
   analysed extinction trials is positive (one-sample t > 0, p < 0.05). The
   change-trial comparison (paper p = 0.053) is reported, not scored.
7. **Background rewards erase the cue response faster than extinction**
   (`fig4/cumsum_bgdrw_vs_extinction.m`). Paired t, n = 7. Pass: same sign as
   paper (t(6) = −3.51), p < 0.05.
8. **Trial-less task: intermediate cue response positive but weaker**
   (`fig5/intermediate_rsp_trialless.m`). Pass: ratio of intermediate to previous
   cue response > 0 (t > 0, p < 0.05; paper t(6) = 6.64), and reward ratio test
   against 1 has the paper's sign (t(6) = 2.95) with p < 0.05.
9. **No backpropagation in trace conditioning** (`fig6/backpropagation_pavlovian.m`).
   Normalised mean (late 2 to 3 s minus early 0 to 1 s) over the first 100 CS+
   trials. Pass: the model group mean lies inside the data's 95% CI.
10. **Sequential conditioning: CS1 and CS2 rise together**
    (`fig6/backpropagation_sequential.py`, control animals, n = 7). Measure:
    mean(CS2 − CS1 AUC 0 to 1 s, first 50 trials) / mean(CS1, trials 150 to 200).
    Pass: the model group mean lies inside the data's 95% CI.
11. **Inhibition from CS2 to reward does not stop CS1 learning**
    (`fig6/learnign_curve_sequential.py`, inhibited animals). Measure: CS1 AUC 0
    to 1.5 s minus AUC −1.5 to 0 s on the last session, over the first-session
    reward response. Pass: positive in every inhibited animal (the paper: "every
    experimental animal ... acquired positive responses to cue1").

Before any model run, the same code computes every test on the real dopamine
data. If our recomputation disagrees with a statistic the paper reports, both
are recorded; the pass rule keeps the paper's reported conclusion as the target.

Secondary (reported, not scored): for every test, whether the model's group
mean lies inside the data's 95% CI.

## Missing sessions

Some days have no NWB file (for example HJ-FP-M2 days 12 to 15 and 44 to 66).
Primary rule: the model gets no experience on those days. Robustness rule:
gaps of up to 5 days inside one protocol block are filled with copies of the
preceding session of that protocol; longer gaps stay empty.

## Search plan

1. Shared parameters: a scrambled Sobol sample of 4,096 points in the box,
   each run on all animals of Tests 1 to 9. Score = number of tests passed, ties
   broken by a continuous loss (for each failed test, distance from its pass
   region in units of the data SD across animals).
2. The 20 best points are refined by Nelder–Mead on the loss with the integer
   choices held fixed.
3. The best 5 shared sets are run on the sequential cohort (Tests 10, 11).
4. If no shared set passes everything, each animal's parameters are refined
   separately from the best shared set, which the claim allows.

## Impossibility side

For each test that no searched model passes, we look for a written argument
that a subfamily of TD-MS cannot pass it (for example: with η = 1 and the laser
term as written, Test 11 requires the integrated inhibition to be smaller than
the reward's discounted value). Such an argument goes in PROOF.md with every
assumption stated. A search that finds nothing is not a proof and is reported
as "not found in the box".

## Robustness checks

- 20 noise seeds per test (primary), and a second block of seeds 100 to 119.
- Missing-session rule: primary versus fill-short-gaps.
- Δ = 0.05 s for the best model.
- Kernel decay 0.3 s and 1.0 s instead of 0.5 s for the best model.

If any of these flips a test for the best model, the status is Contested.

## Addendum A: richer TD classes (frozen 2026-11-25, before any model run)

Added because a negative result for TD-MS alone would leave open that the state
representation was the problem. Same data, same eleven tests, same pass rules,
noise procedure, scale rule, kernel, Δ, inhibition calibration and missing-day
rule as above. Each class keeps one fixed representation and one parameter set
per animal; what is learned inside the representation is learned online from
the same event stream.

**TD-BS: belief-state TD over a learned semi-Markov task model** (after Daw,
Courville and Touretzky 2006; Starkweather et al. 2017). Hidden states: ITI, and
for each cue type c a chain of wait sub-states wait_{c,k}, k = 1..K, each 0.5 s
wide up to 12 s. Every cue onset opens its own wait instance, so overlapping cues
(the trial-less task) superpose. The agent learns, by counts with forgetting
rate κ per occurrence, the probability π_c that cue c is followed by reward and
the histogram F_c of cue-to-reward delays. After τ seconds with no reward, the
belief that instance j of cue c is still waiting is

    p_j(τ) = π_c (1 − F_c(τ)) / (π_c (1 − F_c(τ)) + 1 − π_c),

and the instance's mass p_j(τ) sits on the sub-state whose bin holds τ; a reward
closes the open instance of the matching type with the earliest expected delay.
Features are b_t = [mass on each wait_{c,k} summed over instances, 1 − total
mass clipped to [0, 1] for ITI, optional bias]. Value and learning are the TD-MS
equations with x_t = b_t. Free parameters: α, γ_s, λ, η, w₀, κ (0.01 to 0.5,
log), K fixed at 24.

**TD-SR: TD over learned successor features** (after Dayan 1993; Gardner,
Schoenbaum and Gershman 2018). Base features φ_t are the TD-MS microstimulus
features. A successor matrix M is learned online,

    M ← M + α_M (φ_{t−1} + γ M φ_t − M φ_{t−1}) φ_{t−1}ᵀ,   M initialised to I,

and value is TD(λ) over ψ_t = M φ_t: V_t = w · ψ_t, with the TD-MS rule for δ, the
eligibility trace (over ψ) and the w update. Dopamine is δ. Free parameters: the
TD-MS ones plus α_M (1e−4 to 0.1, log).

Search for each added class: the same procedure as for TD-MS with a 1,024-point
Sobol sample instead of 4,096 (each run costs more), then Nelder–Mead on the 10
best. The decision rule in CLAIM.md applies to the union of the three classes:
one model from any class that passes all eleven tests disproves the statement.
An impossibility argument is stated per class.

## Deviations from plan

None yet.
