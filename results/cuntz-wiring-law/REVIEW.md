# Review

| Date | Reviewer | Person or AI agent | Scope | Outcome |
| :-- | :-- | :-- | :-- | :-- |
| 2026-11-25 | Claude (separate agent session, not the author session) | AI agent | Rebuilt the pipeline from CLAIM.md and METHODS.md before reading the code or outputs; recomputed every β, interval and verdict for both datasets; then compared with the recorded output, reviewed check.py and wiring.py, and checked the git order of the plan and results. | No blocking problems. Same sample gives identical per-neuron values and the same verdicts (FlyWire Refuted, MICrONS Contested, status Contested). Raised: FlyWire's Refuted is marginal and depends on the draw (deviation 3); exclusions applied once at τ = 1 µm (deviation 4); benchmark recipe differs from the text (deviation 5); a vertex left isolated by pruning can enter the hull (V changes ≤ 2.7% in 2.5% of neurons at τ = 5 µm); FlyWire skeletons often have several disconnected pieces (median 2), and restricting to one-piece neurons gives β = 0.535. |

An agent review cannot by itself earn ✅. A person still needs to review this
result.
