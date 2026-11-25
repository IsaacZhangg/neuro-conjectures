# Claim 2.6 at κ = 0.1, 0.25, 0.5, 0.75, 1.0 (certified, Arb)

flint_hess.py (agent-written, transcribes the general-κ Hessian display in the proof of Huang's
Lemma B.9) was run on parameter balls that contain the certified enclosures of (α⋆, q0, ψ0)
from code/tight_*.json (Condition 3.1 certification, α box ±1e-13). Every run prints
"CERTIFIED negative definite on box: True".

| κ | λ_max(∇²S̄⋆(1,0)) | λ_min |
|---|---|---|
| 0.1 | −0.004444085 ± 3e-10 | −0.0555015 |
| 0.25 | −0.004458165 ± 6e-10 | −0.0455896 |
| 0.5 | −0.004344786 ± 4e-10 | −0.0307704 |
| 0.75 | −0.00396175 ± 2e-9 | −0.0185699 |
| 1.0 | −0.00269164 ± 3e-9 | −0.0097232 |

Review items: (1) the transcription of the Lemma B.9 formula (the agent checked it against finite
differences of the float evaluator to 5e-7, far below these margins); (2) the tail constant 1000,
which I checked by hand: every integrand is at most 400(1+|z|)² for q ∈ [0.5, 0.75], ψ ≤ 8, |κ| ≤ 1.
