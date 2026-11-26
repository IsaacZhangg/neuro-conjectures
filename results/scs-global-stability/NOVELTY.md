# Novelty audit

Audit date: 2026-11-25. Done by an AI agent literature search (Firecrawl
paper search with full-text passages, web search; Exa was unavailable), plus
the earlier search recorded in SOURCES.md. Passages were read for the papers
marked "read"; the rest are cited from abstracts or memory, as noted.

## Verdict: partly known

- **For the SCS network itself (tanh, iid Gaussian J), no rigorous global
  stability bound above g = 1/2 was found in print**, at 1/√2 or anywhere up
  to 1. The g < 1 statement is a large-N claim, not a proof.
- **The mechanism is known.** The same two steps (a diagonal Lyapunov
  criterion plus the semicircle edge √2 of the symmetric part) are published
  and explicit for the random Lotka–Volterra model, giving the same constant
  1/√2. A referee could fairly call our Theorem 2 an immediate transfer.
- **What is ours:** stating and proving it for SCS, the non-asymptotic tail
  bound 1 − exp(−N s²/2) valid simultaneously for all g ≤ 1/(√2 + s), the
  explicit convergence rate, and the remark that the identity weight cannot
  go past 1/√2.

So the honest description is: a rigorous push from 1/2 to 1/√2 for SCS by a
known method. The window 1/√2 ≤ g < 1 is the genuinely open part.

## What is known

| Work | What it shows | Rigour | Range |
| :-- | :-- | :-- | :-- |
| Wainrib and Touboul, PRL 110, 118101 (2013), arXiv:1210.5082 (read) | Unique, globally stable equilibrium below σ = 1, via the bulk spectral disc of J·diag(S′(x)) and a contraction argument | Heuristic: the bulk disc does not bound every eigenvalue uniformly in x, and Hurwitz Jacobians do not by themselves give global stability (Markus–Yamabe) | claims g < 1 |
| Sherf, Tang, Hafner, Touboul, Pitkow, Bassler, Josić, arXiv:2512.24439 (read) | Partially symmetric SCS; says the system is contracting for g(1+τ) < 1 | Argument not given; uniqueness proof assumes hyperbolic stable fixed points; spectral bound is large-N | claims g(1+τ) < 1 |
| Stubenrauch et al., Phys. Rev. Research 7, 023203 (2025), arXiv:2210.07877 (read) | Fixed-point counting for g > 1 | Kac–Rice plus replica, physics level; cites Wainrib–Touboul for g < 1 | g > 1 only |
| Kadmon and Sompolinsky, PRX 2015; Crisanti and Sompolinsky, PRE 2018 | Transition to chaos at g = 1 | Dynamic mean-field theory, large N | g ≈ 1 |
| Touboul, arXiv:2610.38893 (review) | Rigorous mean-field results for such networks hold on fixed time windows; uniform-in-time questions open | Survey | n/a |
| Ben Arous and Guionnet 1995; Dembo, Lubetzky and Zeitouni, arXiv:1911.08001 | Rigorous large-N limits of the dynamics | Rigorous, finite time horizon only | n/a |
| Matsuoka 1992; Forti and Tesi 1995 (from memory) | Lyapunov diagonal stability of J − I implies global stability for x′ = −x + Jφ(x) | Rigorous, deterministic | any J |
| Kozachkov et al., PLoS CB 2020 and arXiv:2106.08928; Centorrino et al., arXiv:2302.13452; Davydov, Proskurnikov and Bullo, arXiv:2110.08298 (read) | Contraction criteria for rate networks (symmetric-part, LDS-type, ℓ1/ℓ∞) | Rigorous, deterministic; not applied to Gaussian J, no g threshold stated | any J |
| Clenet, Massol and Najim, J. Math. Biol. 2023, arXiv:2205.00735, Thm 2 (read) | Random Lotka–Volterra: unique globally stable equilibrium when the interaction strength is below 1/√2 (their α > √2), via Takeuchi–Adachi plus the semicircle edge | Rigorous | same constant 1/√2, different model |
| Clenet, El Ferchichi and Najim, arXiv:2205.15591, Prop. 2.3 (read) | Elliptic version: σ < 1/√(2(1+ρ)) | Rigorous | different model |
| Gueddari, Hachem and Najim, arXiv:2402.08271, Remark 3.5 (read) | States that beyond the √(2(1+ρ)) condition existence of a stable equilibrium "is not granted" | Says the next step is open | different model |
| Fyodorov and Khoruzhenko, PNAS 2016; Ben Arous, Fyodorov and Khoruzhenko, PNAS 2021 | Equilibrium counting for related random fields | Expected counts; do not give global convergence | different model |
| Cessac, J. Phys. A 27, L927 (1994) (abstract only) | Absolute-stability gain threshold for a discrete-time random network | Large-N | likely the operator-norm analogue |

## What is open

- Global stability with high probability for 1/√2 ≤ g < 1 in the SCS model.
  This needs a non-identity weight (Lyapunov diagonal stability of
  (g/√N)G − I for Ginibre G past 1/√2, which we did not find studied), a
  non-quadratic Lyapunov function, or an argument that uses the dynamics.
- The same step past 1/√2 is open for random Lotka–Volterra (Gueddari,
  Hachem and Najim, Remark 3.5).

## Effect on the catalogue

The result stays Partial for the conjecture (proved for g < 1/√2, open for
1/√2 ≤ g < 1). README.md and the catalogue row now say that the method is
known from the Lotka–Volterra literature.

## Addendum 2026-11-26: the diagonal-weight barrier (BARRIER.md)

**Verdict: new as far as we can find.** BARRIER.md, Theorem 3, proves that for
every fixed g > 1/√2, J − I = (g/√N)G − I is not Lyapunov diagonally stable
with probability ≥ 1 − exp(−cN), so the threshold for diagonal stability of
Ginibre-type matrices is exactly 1/√2 (and 1/√(2(1+ρ)) for zero-mean Gaussian
elliptic pairs). The question was posed as open by Gibbs, Grilli, Rogers and
Allesina, Phys. Rev. E 98, 022410 (2018), arXiv:1708.08837, Section VII
(read; both quoted sentences checked in the paper): "Under which conditions
they are also diagonally stable is an important open problem." A second agent
search (Firecrawl paper index and web, 2014–2026; Exa unavailable) found no
paper that proves or refutes it. Searched and not relevant: Clenet–Massol–Najim
arXiv:2205.00735 and Akjouj–Hachem–Maïda–Najim arXiv:2302.07820 (use D = I,
do not discuss optimality); Akjouj et al. arXiv:2212.06136 (say the √2
condition "might not be necessary" without addressing diagonal weights);
Chen–Wang–Liu, Phys. Rep. 2024, arXiv:2312.07737; Stone 2018 (D-stability
only); Clenet–Gueddari arXiv:2603.01774; Hachem arXiv:2302.09847;
Gueddari–Hachem arXiv:2510.15754; Bizeul–Najim arXiv:1904.04559;
Ben Arous–Fyodorov–Khoruzhenko arXiv:2008.00690; Kozachkov and
Davydov–Bullo contraction papers. The proof's ingredients (the easy direction
of Barker–Berman–Plemmons 1978, GOE invariance, independence of symmetric and
antisymmetric Gaussian parts) are classical; the combination and the result
appear to be new. This does not change the conjecture's status: the window
1/√2 ≤ g < 1 is still open, and Theorem 3 shows that closing it needs a
Lyapunov function outside the diagonal family.

