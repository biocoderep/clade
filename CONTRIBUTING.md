# Contributing to CLADE

## Development setup

```
git clone https://github.com/biocoderep/clade.git
cd clade
python -m venv .venv && source .venv/bin/activate
pip install -e ".[firth,dev]"
pytest tests/
```

## Before opening a PR

- `pytest tests/` must pass.
- `ruff check src/ tests/` must pass.
- If you're changing the statistical logic in any of the six validation stages (`src/clade/validation/`), explain why in the PR description — these are meant to be stable, auditable implementations, not code that changes casually.
- New dependencies go in `pyproject.toml`, `environment.yml`, and `requirements.txt` together — a dependency declared in only one of these was a real bug caught during this project's own development (see `docs/reproducibility/open_items.md`).

## Reporting a candidate that CLADE misclassifies

Open an issue with the `bug_report` template, including the exact `StageResults` values (or the CSV inputs) that produced the unexpected disposition. `src/clade/classification/disposition.py`'s docstring is explicit that the rule-based classifier codifies the *general* case, not every nuanced judgment call from the original manual evidence table — genuinely borderline candidates are expected to need a second look, not treated as bugs by default.

## Extending CLADE to a new organism or resistance mechanism

This has not been done yet — CLADE has been validated on one organism, one resistance determinant, and one dataset (see `CLADE_Framework_Specification.md`, Section 6). If you try it on a new dataset, please open an issue either way (success or failure) — that's exactly the kind of report this project needs to become more than a single case study.
