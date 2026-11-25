# Claim

## Statement tested

Random rate network of Sompolinsky, Crisanti and Sommers (SCS):

    dx/dt = -x + J tanh(x),   x in R^N,   J_ij iid N(0, g^2/N), no other structure
    (self-couplings J_ii included, drawn like the rest).

**Conjecture (part a of target #8).** For every fixed g < 1, the origin is
globally asymptotically stable (every solution converges to 0, and 0 is
Lyapunov stable) with probability tending to 1 as N tends to infinity.

This is our precise reading of the folklore statement "for g < 1 the network
decays to the silent fixed point". The physics literature derives it from
dynamical mean-field theory, which only describes typical initial conditions
and finite time windows; the conjecture asks for all initial conditions.

## Where it comes from

- Sompolinsky, Crisanti and Sommers, "Chaos in random neural networks",
  Phys. Rev. Lett. 61, 259 (1988): the zero fixed point is the stable state
  for g < 1 and chaos appears for g > 1.
- Target #8 in the conjecture-finder candidate list (2026-11-25): "the spectral
  norm of J is about 2g, which only gives a contraction proof for g < 1/2".

## What counts

This is a theorem target. Numerics guide the work and never count as the
result.

- **Proved** if a complete written proof covers every g < 1.
- **Partial** if a complete written proof covers g < g_1 for some explicit
  g_1 > 1/2 but not all g < 1. The README must state g_1 and the
  probability bound.
- **Open** if no complete proof beyond g < 1/2 is obtained.

Decision rule frozen in this commit, before the proof work is written up.
