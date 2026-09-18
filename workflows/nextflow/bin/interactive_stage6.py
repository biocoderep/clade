#!/usr/bin/env python3
"""Interactive Stage 6: external/literature corroboration.

Stage 6 is deliberately NOT automated anywhere in CLADE. It is a judgement
about mechanistic plausibility -- reading the literature, checking homology,
deciding whether a proposed function makes biological sense -- and a pipeline
that pretended to automate that would be asserting evidence it had not
gathered. Every other stage in this framework runs unattended inside the
Nextflow DAG; this is the one stage that does not, on principle, and this
script is where that human judgement actually happens.

Run it after the core evidence table exists (RUN_CLADE has produced
clade_evidence_table.tsv). It shows you each candidate's disposition so far --
what to prioritise is obvious from that, not guessed at -- and asks, for each
one, whether you found external mechanistic or literature support. Answers
are appended to (or start) a stage6.json you then pass back into a second
pipeline run via --stage6, which only touches Stage 6 bookkeeping: it cannot
change a REJECTED or INSUFFICIENT_EVIDENCE verdict, because a disposition that
already required an explicit pass at Stages 1/3/4/5 does not become
convergent on the strength of a literature search alone (CLADE's own rule,
enforced by clade.classification.disposition, not by this script).

Usage:
    interactive_stage6.py --evidence-tsv clade_evidence_table.tsv \\
        --output stage6.json [--only-unresolved] [--resume]

Answers accepted at each prompt: y / n / u (unknown -- skip, ask again later)
"""
from __future__ import annotations

import argparse
import json
import sys

import pandas as pd


def prompt_one(candidate: str, row: pd.Series) -> bool | None:
    print(f"\n{'=' * 70}")
    print(f"Candidate: {candidate}")
    print(f"  Disposition so far : {row.get('Disposition', '?')}")
    stage_cols = [c for c in row.index if c.endswith("_status")]
    for c in stage_cols:
        print(f"  {c:<16}: {row[c]}")
    reason = row.get("Reason")
    if isinstance(reason, str) and reason.strip():
        print(f"  Reason             : {reason}")
    print(f"{'=' * 70}")
    while True:
        ans = input(
            "Did you find external mechanistic/literature support for this "
            "candidate?  [y/n/u=unknown, s=skip, q=quit] "
        ).strip().lower()
        if ans in ("y", "yes"):
            return True
        if ans in ("n", "no"):
            return False
        if ans in ("u", "unknown", ""):
            return None
        if ans in ("s", "skip"):
            return "__skip__"
        if ans in ("q", "quit"):
            return "__quit__"
        print("  Please answer y, n, u, s, or q.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidence-tsv", required=True,
                    help="clade_evidence_table.tsv from a completed RUN_CLADE pass")
    ap.add_argument("--output", required=True, help="stage6.json to write")
    ap.add_argument("--only-unresolved", action="store_true",
                    help="prompt only for candidates not already REJECTED or UNINFORMATIVE "
                         "-- Stage 6 cannot change either of those verdicts, so most users "
                         "want this on")
    ap.add_argument("--resume", action="store_true",
                    help="load --output if it already exists and skip candidates already "
                         "answered in it")
    args = ap.parse_args()

    if not sys.stdin.isatty():
        sys.exit(
            "interactive_stage6.py needs a real terminal (stdin is not a tty here).\n"
            "Run it directly in your shell, not inside a Nextflow process or a "
            "non-interactive script -- Stage 6 is deliberately the one stage in "
            "this framework that is not automated."
        )

    table = pd.read_csv(args.evidence_tsv, sep="\t").set_index("Candidate")

    answers: dict[str, bool | None] = {}
    if args.resume:
        try:
            with open(args.output) as fh:
                answers = json.load(fh)
            print(f"Resuming: {len(answers)} candidate(s) already answered in {args.output}")
        except FileNotFoundError:
            pass

    candidates = list(table.index)
    if args.only_unresolved:
        skip_dispositions = {"rejected", "uninformative"}
        candidates = [
            c for c in candidates
            if str(table.loc[c, "Disposition"]).lower() not in skip_dispositions
        ]
        print(f"--only-unresolved: {len(candidates)} of {len(table)} candidates need Stage 6 "
              f"(REJECTED and UNINFORMATIVE candidates cannot become convergent on Stage 6 "
              f"alone, and are skipped).")

    todo = [c for c in candidates if c not in answers]
    if not todo:
        print("Nothing left to answer.")
    else:
        print(f"\n{len(todo)} candidate(s) to review. Answers save after every candidate, "
              f"so it is safe to quit (q) and resume later with --resume.\n")

    for c in todo:
        row = table.loc[c]
        result = prompt_one(c, row)
        if result == "__quit__":
            break
        if result == "__skip__":
            continue
        answers[c] = result
        with open(args.output, "w") as fh:
            json.dump(answers, fh, indent=2, sort_keys=True)

    with open(args.output, "w") as fh:
        json.dump(answers, fh, indent=2, sort_keys=True)

    n_yes = sum(1 for v in answers.values() if v is True)
    n_no = sum(1 for v in answers.values() if v is False)
    n_unk = sum(1 for v in answers.values() if v is None)
    print(f"\nWrote {args.output}: {len(answers)} recorded "
          f"({n_yes} support found, {n_no} searched and found none, {n_unk} unknown).")
    print("Re-run the pipeline with --stage6 " + args.output + " to fold these into "
          "the evidence table. This changes Stage 6 bookkeeping only -- it cannot promote "
          "a REJECTED or INSUFFICIENT_EVIDENCE candidate to CONVERGENT.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
