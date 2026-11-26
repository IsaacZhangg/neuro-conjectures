# Paused 2026-11-26 (screening, 2 noise seeds)

Work was paused on request. Saved here so a later session can resume.

| File | Points | State |
|---|---|---|
| search_ms.jsonl.gz | 2,957 of 4,096 | TD-MS Sobol screen, Tests 1 to 9; best 7 of 9 (fails 7, 9) |
| search_bs.jsonl.gz | 1,024 of 1,024 | TD-BS screen done; best 6 of 9; Test 1 never passes |
| search_sr.jsonl.gz | 905 of 1,024 | TD-SR screen; best 6 of 9 (fails 7, 8, 9) |
| refine_bs.jsonl.gz | 10 starts x 60 | Nelder-Mead on the 10 best TD-BS points; best 7 of 9 (fails 1, 5) |
| seq_bs.jsonl.gz | 2 | Two best refined TD-BS sets on all cohorts: 8 of 11 (fail 1, 5, 10) |

These are screening numbers with 2 noise seeds, not results.

To resume: gunzip into the scratch folder and rerun `search.py` with the same
`--n`, `--seed` and `--cls`; it skips indices already in the output file.
Then `refine.py --cls MS --top 20` and `--cls SR --top 10`, then
`seqcheck.py` on the best five across classes, then the 20-seed rescoring and
robustness checks in METHODS.md.
