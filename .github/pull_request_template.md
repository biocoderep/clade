## What this changes

## Why

## Checklist

- [ ] `pytest tests/` passes locally
- [ ] `ruff check src/ tests/` passes locally
- [ ] If this changes a validation stage's logic, I've explained why in the PR description — Stages 1–5's statistical methods should not change silently (see `CLADE_Framework_Specification.md`)
- [ ] If this adds a new dependency, it's declared in `pyproject.toml` (this is exactly the kind of gap that caused a real CI-catchable bug during development — see `docs/reproducibility/open_items.md`)
- [ ] No real dataset files were added (see `.gitignore` and `DATA_AVAILABILITY.md`)
