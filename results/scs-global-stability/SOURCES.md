# Sources

Checked 2026-11-25.

## Primary

- H. Sompolinsky, A. Crisanti, H. J. Sommers. Chaos in random neural networks.
  Phys. Rev. Lett. 61, 259 (1988). DOI 10.1103/PhysRevLett.61.259. The silent
  fixed point is the stable state for g < 1; chaos for g > 1 (mean-field
  theory, not a proof of global stability).

## Results the proof uses

- K. Matsuoka. Stability conditions for nonlinear continuous neural networks
  with asymmetric connection weights. Neural Networks 5, 495–500 (1992).
  DOI 10.1016/0893-6080(92)90011-7. Lyapunov-diagonal-stability criterion for
  x' = −x + Wf(x); Theorem 1 here is its case D = I.
- M. Forti, A. Tesi. New conditions for global stability of neural networks
  with application to linear and quadratic programming problems. IEEE Trans.
  Circuits Syst. I 42, 354–366 (1995). Same criterion, class of sigmoids.
- Sudakov–Fernique inequality: R. J. Adler, J. E. Taylor, Random Fields and
  Geometry (Springer 2007), Thm 2.2.3.
- Gaussian concentration for Lipschitz functions: S. Boucheron, G. Lugosi,
  P. Massart, Concentration Inequalities (OUP 2013), Thm 5.6.
- Edge of Wigner matrices: Z. D. Bai, Y. Q. Yin, Ann. Probab. 16, 1729 (1988);
  Z. Füredi, J. Komlós, Combinatorica 1, 233 (1981).

## Later work that retests, refines or bears on it

- H. Suetani. Broad onset of chaos with toroidal dynamics in finite-size
  random neural networks. Chaos 36, 093107 (2026). DOI 10.1063/5.0341330.
  Numerical; reports no non-origin attractors below g = 1 at finite N; no
  proof of global stability.
- D. Pazó. Discontinuous transition to chaos in a canonical random neural
  network. Phys. Rev. E 110, 014201 (2024), arXiv:2405.14607. With
  φ = tanh + ε tanh³, chaos can appear below g = 1. For those ε the slope of φ
  exceeds 1, so it does not contradict Theorem 1, but it shows that "g < 1
  implies global stability" depends on the nonlinearity.
- G. Wainrib, J. Touboul. Topological and dynamical complexity of random
  neural networks. Phys. Rev. Lett. 110, 118101 (2013). Mean number of
  equilibria jumps from one to exponentially many at g = 1 (not rigorous).
- J. Touboul. Universality in randomly connected neural networks: old and new
  on a conjecture of Geman and Hwang. arXiv:2610.38893 (2026). States that all
  rigorous mean-field results hold on fixed time intervals only.
- A. Davydov, A. Proskurnikov, F. Bullo. Non-Euclidean contraction analysis of
  continuous-time neural networks. arXiv:2110.08298. Weighted l1/l∞ tests;
  for iid Gaussian J these grow with √N and give nothing better (our reading).
- Gokhale, Proskurnikov, Kawano, Bullo. arXiv:2604.00119 (2026). Sharp
  contraction LMIs for firing-rate networks; Lyapunov diagonal stability of
  W − I is necessary for their test, so contraction certificates do not beat
  the diagonal threshold (no random-matrix analysis there).
- Random Lotka–Volterra analogue (same gap between the symmetric-part
  condition and the expected threshold): M. Clenet, H. El Ferchichi,
  J. Najim, arXiv:2205.15591; M. Clenet, F. Massol, J. Najim, J. Math. Biol.
  (2023), arXiv:2205.00735; I. Akjouj et al., review, arXiv:2212.06136.

## Searches run

- Firecrawl paper index: "global asymptotic stability random recurrent neural
  network tanh Gaussian coupling g<1 proof Lyapunov"; "random recurrent neural
  network x' = -x + J tanh(x) fixed point globally stable for coupling below
  one, rigorous" (from 2015).
- Web: "Sompolinsky Crisanti Sommers random network global stability origin
  g<1 rigorous proof"; "Lyapunov diagonal stability random matrix sqrt(2)
  threshold Lotka-Volterra global stability unique equilibrium".
- An AI research helper also searched contraction-analysis papers
  (Kozachkov–Slotine, Miller–Hardt, Centorrino et al., Davydov–Bullo),
  rigorous DMFT (Ben Arous–Guionnet, Dembo et al., Celentano–Cheng–Montanari),
  and Zames–Falb / IQC work on repeated nonlinearities. Some sources were read
  only as abstracts or excerpts.
- We found no paper that states global stability for SCS beyond g < 1/2, and
  none that states the g < 1/√2 corollary.
