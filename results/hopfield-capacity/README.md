# Hopfield network storage capacity

A Hebbian Hopfield network can recall stored random patterns only while the
number of patterns stays below about 0.138 times the number of neurons.

Status: Supported
Kind: Model claim (replication)
Reproducibility: 🟡 runnable check with recorded output; not yet independently reviewed

## Result

Simulated networks of N = 500, 1000, 2000 and 4000 neurons. The load at which
retrieval succeeds half the time falls with N (0.172, 0.163, 0.155, 0.153).
Extrapolated to large N with an N^-1/2 correction, the capacity is **0.141**,
within the pre-set tolerance of 0.138 ± 0.02. Below the transition the final
overlap matches the theory's error rates closely (for example, mean overlap
0.998 at α = 0.10, where the theory predicts about 0.1% errors).

## Files

- [CLAIM.md](CLAIM.md): the exact statement and decision rule.
- [SOURCES.md](SOURCES.md): primary references.
- [METHODS.md](METHODS.md): model, dynamics, parameters, limits.
- [check.py](check.py): the simulation. `python check.py` (about 4 minutes on
  a multi-core laptop) or `python check.py --quick` (under a second).
- [check-output.json](check-output.json): recorded full run.
- [REVIEW.md](REVIEW.md): review record (pending).

## Scope and limits

- Tests the model, not the brain. Real cortex is not a symmetric Hebbian network
  of binary units.
- Finite sizes only. The large-N value depends on an assumed N^-1/2 correction;
  the fit uses 4 sizes with 120 samples per point, so the N = 2000 and N = 4000
  estimates are within sampling noise of each other.
- Tests replica-symmetric capacity at zero temperature with retrieval defined as
  final overlap ≥ 0.9. Other definitions shift the finite-size numbers.
- This is a replication of a well-established result. It is included as the
  worked example for the collection's format, not as a new finding.
