# Literature notes: global stability of x' = -x + J tanh(x), J = (g/sqrt N) G, routes past g = 1/sqrt2
(compiled 2026-11-26; "READ" = read full-text passages; "ABS" = abstract/metadata only; "KNOWN" = from background knowledge, not re-verified this session)

## Q2. Random rate networks: rigorous results

### Touboul, "Universality in randomly connected neural networks: old and new on a conjecture of Geman and Hwang", arXiv:2610.38893 (Oct 2026). READ (passages)
- Proves: (weakly heterogeneous, Var o(1/n)) quenched convergence to McKean-Vlasov on every finite time interval under 4th moment; (strongly heterogeneous, Var ~ 1/n, our regime) convergence of the network to the DMFT/self-consistent Gaussian-process limit, AVERAGED over sub-Gaussian couplings, on every FINITE time interval [0,T], via relative-entropy comparison to Gaussian couplings; heavy-tailed alpha-stable conjecture, existence/uniqueness of limit.
- Nothing on long-time behaviour / global stability / convergence to 0 for g<1. All results are finite-horizon. It explicitly notes rare large couplings "may affect long-time stability and the transition to chaos" (heavy tails) — i.e., long-time is outside scope.
- Useful review pointers: Geman-Hwang (linear S, universality), Ben Arous-Guionnet 1995/1997, Cabana-Touboul 2013/2018 (LDP for random neural nets with Gaussian couplings, finite T), Celentano-Cheng-Montanari arXiv:2112.07572.

### Sherf, Tang, Hafner, Touboul, Pitkow, Bassler, Josic, "Complexity and dynamics of partially symmetric random neural networks", arXiv:2512.24439 (Dec 2025). READ (appendices A, B)
- Model x' = -x + W S(x), elliptic ensemble with reciprocity tau, g_eff = g(1+tau).
- Claims "when g_eff<1 the system is contracting as N->inf and x*=0 is the unique fixed point". The argument is NOT a rigorous global-stability proof:
  * App. A: rightmost eigenvalue of W*Delta(d) <= g(1+tau) is derived from the deterministic-equivalent / inhomogeneous-elliptic-law equations of ref [21] (Kuczala-Sharpee-type), which hold for a FIXED profile d independent of W, not uniformly over d chosen adversarially (d = S'(x*) depends on W).
  * App. B: Brouwer degree argument: IF every fixed point has Hurwitz Jacobian, then the fixed point is unique (degree = (-1)^N). This is uniqueness of equilibria, not global attractivity (does not exclude limit cycles/chaos).
- Takeaway: at most heuristic evidence for g<1 uniqueness; the degree trick is a reusable lemma: "Hurwitz at every equilibrium + homotopy to (g_eff-1)x  =>  unique equilibrium".
