# Neuro Conjectures

Reproducible tests of conjectures, predictions and hypotheses in neuroscience.
Each result states one exact claim from the literature, tests it with code that
anyone can rerun, and labels what the evidence does and does not show.

Neuroscience claims are mostly empirical, so "verified" here does not mean
"proved". A result is checked in one of three ways:

- **Theorem**: a mathematical statement from theoretical neuroscience, checked
  by a written proof and, where practical, a Lean 4 formalization.
- **Model claim**: a statement about a model (a network, a cable equation, a
  coding scheme), checked by deterministic simulation with recorded output.
- **Empirical claim**: a statement about real nervous systems, checked by a
  pre-registered reanalysis of public data (Allen Brain Observatory, DANDI,
  OpenNeuro, NeuroMorpho.org, CRCNS and similar).

Curated by Isaac Zhang with Claude (Anthropic) as AI co-curator. External
contributions are credited on each result.

All results are unrefereed. Internal review is not peer review. Nothing here is
medical advice.

[Catalogue](CONTENTS.md) · [How results are checked](docs/verification.md) ·
[Collection design](DESIGN.md) · [Candidate targets](docs/candidates.md) ·
[History](HISTORY.md)

## Results

| # | Claim | Kind | Status | Reproducibility |
| :-- | :-- | :-- | :-- | :-: |
| 1 | [Hopfield network storage capacity, α_c ≈ 0.138](results/hopfield-capacity/) | Model claim (replication) | Supported, finite-size simulation | 🟡 |
| 2 | [Neuropeptide links and C. elegans signal propagation](results/celegans-neuropeptides/) | Empirical claim | Partial: helps only for pairs with no wired path, one peptide model | 🟡 |
| 3 | [All-pattern capacity of dense associative memory, Gaussian patterns](results/dense-am-capacity/) | Theorem | Proved (static retrieval criterion): α_c = union bound for all λ, (log 2)/4 for λ ≥ 0.70092 | 🟡 |
| 4 | [SCS random network: silent state globally stable for g < 1/√2](results/scs-global-stability/) | Theorem | Partial: proved for g < 1/√2 (method known from random Lotka–Volterra), open for 1/√2 ≤ g < 1 | 🟡 |
| 5 | [Rigorous Hopfield energy-barrier bound, α ≥ 0.1295 for ±1 patterns](results/hopfield-rigorous-bound/) | Theorem | Partial: ±1 barrier proof written, passed independent agent review (previous rigorous barrier bound 0.08); Gaussian cases proved; 0.1347 fixed point superseded if Zhou (2026, unrefereed) holds | 🟡 |
| 6 | [Tree committee machine, K = 3: capacity at least 1.55 (new proved lower bound)](results/tree-committee-k3/) | Theorem | Partial, lower bound only, Gaussian patterns | 🟡 |

Status words are defined in [DESIGN.md](DESIGN.md#status). Reproducibility marks:

- ✅ rerun by CI from a clean environment **and** independently reviewed
- 🟢 rerun by CI from a clean environment
- 🟡 runnable check with recorded output, not yet in CI or not yet reviewed
- ⬜ written argument or analysis plan only

## Layout

```
results/<slug>/      one folder per claim (see templates/result/)
templates/result/    the files every result starts from
docs/verification.md what each kind of check must contain
docs/candidates.md   screened targets not yet attempted
scripts/             repository checks (catalogue consistency)
```

## Contributing

Open an issue to claim a target before starting, using the screening questions
in [DESIGN.md](DESIGN.md#screening-a-target). Copy `templates/result/` to
`results/<slug>/`, freeze your analysis plan in a commit before you touch the
data, then add the check and its recorded output. Run
`python scripts/check_catalogue.py` before opening a pull request.

## How to cite

See [CITATION.cff](CITATION.cff). For a specific result, also cite the original
paper that made the claim, which is listed in that result's `SOURCES.md`.

## License

Code is under the MIT License ([LICENSE](LICENSE)). Text, figures and recorded
outputs are under CC BY 4.0. Each dataset keeps its own license, recorded in the
result that uses it.
