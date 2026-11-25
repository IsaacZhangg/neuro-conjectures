# Instructions for coding agents

These rules apply to any AI agent working in this repository.

- Read DESIGN.md and docs/verification.md before starting.
- Claim a target in an issue before working on it, and check the catalogue and
  open issues for the same claim under another name.
- Never change a claim's wording, quantities or conditions to make it pass.
  If the original claim is vague, write the precise version in CLAIM.md and say
  that it is your reading.
- For empirical claims, commit METHODS.md and the decision rule before any code
  reads the data. Do not tune the analysis after seeing results without
  recording it under "Deviations from plan".
- Record what you could not check. "Not verified" is an acceptable line;
  an unsupported "verified" is not.
- Do not commit raw data, credentials, or anything that could identify a person.
- An agent review is labelled as such in REVIEW.md and cannot by itself earn ✅.
- Before a pull request: run every touched check in quick mode
  (`scripts/quick_checks.sh` runs them the way CI does, one uv environment per
  result built from its own requirements.txt) and
  `python scripts/check_catalogue.py`, and update README.md, CONTENTS.md and
  HISTORY.md together.
