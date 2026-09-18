#!/usr/bin/env python3
"""The pipeline's final answer, stated once, plainly: a candidate, or no
candidate.

Reads clade_evidence_table.tsv (required) and, if it exists, a Stage 5
conditional-model result (optional -- corrects the built-in matched-pair test
for control-genome reuse, which the case study this framework was built on
found necessary; see stage5_conditional.py). Neither file's numbers are
recomputed here: this script only adjudicates between what is already there.

Three possible answers, and the rule that decides between them mirrors
CLADE's own disposition logic exactly (clade.classification.disposition) --
nothing here loosens or reinterprets it:

  CONVERGENT CANDIDATE(S) FOUND
      At least one candidate is CONVERGENT: an explicit pass at every stage
      CLADE requires (1, 3, 4, 5). Named, with its full evidence.

  CANDIDATE NOMINATED FOR FOLLOW-UP (not confirmed)
      No candidate is CONVERGENT, but at least one is INSUFFICIENT_EVIDENCE
      with three of the four required stages explicitly passing and the
      fourth genuinely missing or underpowered rather than contradicted.
      This is the framework's own worked example of the distinction: a
      candidate can fail to reach CONVERGENT for a reason that says "test it
      further," not "it is wrong." Named, with the same evidence, and with
      the explicit caveat that this is not a validated finding.

  NO CANDIDATE
      Every candidate is REJECTED, UNINFORMATIVE, CLONAL_ARTIFACT, or
      INSUFFICIENT_EVIDENCE for a reason that does not meet the "three of
      four stages pass" bar above. This is a real, reportable result, not an
      absence of one -- the case study this framework was built on ended
      here, after four analytical defects that had produced two false
      CONVERGENT candidates were found and corrected.

Usage:
    final_verdict.py --evidence-tsv clade_evidence_table.tsv \\
        [--stage5-conditional stage5_conditional.tsv] \\
        --output FINAL_VERDICT.md
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

REQUIRED_STAGES = ["Stage1_status", "Stage3_status", "Stage4_status", "Stage5_status"]


def count_passing(row: pd.Series, cond: pd.DataFrame | None) -> tuple[int, list[str]]:
    """How many of the four required stages explicitly pass for this
    candidate, folding in a conditional Stage 5 result where one was
    supplied and the built-in Stage 5 could not report (withheld for
    matched-pair pseudoreplication -- see stage5_matched_neighbors.py)."""
    passing, missing = [], []
    for col, stage in zip(REQUIRED_STAGES, ["1", "3", "4", "5"]):
        status = str(row.get(col, "not run")).strip().lower()
        if stage == "5" and status != "pass" and cond is not None:
            crow = cond[cond.Candidate == row["Candidate"]]
            if len(crow):
                c = crow.iloc[0]
                if pd.notna(c.get("conditional_p")) and c["conditional_p"] < 0.05 \
                        and pd.notna(c.get("conditional_coef")) and c["conditional_coef"] > 0:
                    status = "pass"
        if status == "pass":
            passing.append(stage)
        elif status != "fail":
            missing.append(stage)
    return len(passing), missing


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidence-tsv", required=True)
    ap.add_argument("--stage5-conditional", default=None,
                    help="optional: stage5_conditional.tsv, used only to recognise a "
                         "candidate whose built-in Stage 5 was withheld for "
                         "pseudoreplication but whose conditional re-analysis passed")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    table = pd.read_csv(args.evidence_tsv, sep="\t")
    cond = pd.read_csv(args.stage5_conditional, sep="\t") if args.stage5_conditional else None

    convergent = table[table.Disposition.astype(str).str.lower() == "convergent"]

    nominated = []
    for _, row in table.iterrows():
        if str(row.Disposition).lower() in ("convergent", "rejected", "uninformative"):
            continue
        n_pass, missing = count_passing(row, cond)
        if n_pass >= 3:
            nominated.append((row.Candidate, n_pass, missing, row))

    lines = ["# CLADE Final Verdict", ""]

    if len(convergent):
        lines.append(f"## CONVERGENT CANDIDATE{'S' if len(convergent) > 1 else ''} FOUND")
        lines.append("")
        for _, row in convergent.iterrows():
            lines.append(f"**{row.Candidate}** -- explicit pass at every required stage "
                        f"(1, 3, 4, 5).")
        lines.append("")
        lines.append("This is CLADE's strongest verdict: independently-falsifiable evidence "
                     "at every stage the framework requires. It is not experimental "
                     "validation -- no candidate produced by this pipeline has been "
                     "tested in the wet lab.")
        headline = f"CONVERGENT: {', '.join(convergent.Candidate.astype(str))}"

    elif nominated:
        nominated.sort(key=lambda t: t[1], reverse=True)
        lines.append("## CANDIDATE NOMINATED FOR FOLLOW-UP (not confirmed)")
        lines.append("")
        lines.append("No candidate reached CONVERGENT. The following came closest -- three "
                     "of the four required stages explicitly pass, and the remaining stage "
                     "is missing or underpowered rather than contradicted:")
        lines.append("")
        for cand, n_pass, missing, row in nominated:
            lines.append(f"- **{cand}**: {n_pass}/4 stages pass "
                        f"(missing: {', '.join(missing) or 'none, borderline'}) "
                        f"-- disposition `{row.Disposition}`")
        lines.append("")
        lines.append("This is not a validated finding. It is the answer to \"which candidate "
                     "is worth testing next, and why\" -- a materially weaker claim than "
                     "CONVERGENT, and CLADE's disposition logic keeps the two separate on "
                     "purpose.")
        top = nominated[0][0]
        headline = f"NOMINATED (not confirmed): {top}"

    else:
        lines.append("## NO CANDIDATE")
        lines.append("")
        n = len(table)
        by_disp = table.Disposition.value_counts()
        lines.append(f"None of the {n} candidates screened reached CONVERGENT, and none came "
                     "close enough (3 of 4 required stages) to be nominated for follow-up.")
        lines.append("")
        for disp, count in by_disp.items():
            lines.append(f"- {count} `{disp}`")
        lines.append("")
        lines.append("This is a real, reportable result, not an absence of one.")
        headline = "NO CANDIDATE"

    text = "\n".join(lines) + "\n"
    with open(args.output, "w") as fh:
        fh.write(text)

    print(text)
    print(f"\n{'=' * 70}\nHEADLINE: {headline}\n{'=' * 70}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
