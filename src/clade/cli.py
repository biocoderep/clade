"""
CLADE command-line interface.

    clade --version
    clade --help
    clade validate --genotypes G1.csv [G2.csv ...] --phenotype pheno.tsv \\
        --lineages mlst.tsv --candidates cand1 cand2 ... \\
        [--tree tree.nwk] [--distances distances.tsv] \\
        --output evidence_table.md

Stages 4 (temporal ordering) and 5 (matched-neighbor) are skipped, with a
clear message, if `--tree` / `--distances` are not supplied — CLADE does
not fabricate a result for a stage it wasn't given the input to run.
Stage 6 (external corroboration) is never run by the CLI: it is a manual
literature/homology process by design (see
clade.validation.stage6_corroboration), not something to automate here.
"""
from __future__ import annotations

import argparse
import sys

from clade import __version__


def _cmd_validate(args: argparse.Namespace) -> int:
    import dendropy
    import pandas as pd

    from clade.classification.disposition import StageResults
    from clade.io.genotypes import load_joined_dataset, st_covariate_dummies
    from clade.reporting.evidence_table import build_evidence_table, evidence_table_to_markdown
    from clade.validation.stage1_association import firth_association, fisher_fdr_screen
    from clade.validation.stage2_recurrence import lineage_recurrence
    from clade.validation.stage3_direction import direction_check
    from clade.validation.stage4_temporal_order import cohort_baseline_prevalence, temporal_ordering
    from clade.validation.stage5_matched_neighbors import matched_neighbor_test, nearest_negative_neighbor

    df = load_joined_dataset(args.phenotype, args.genotypes, args.lineages)
    pheno_col = "phenotype"
    candidates = args.candidates

    print(f"[CLADE] Loaded {len(df)} samples, {len(candidates)} candidates.", file=sys.stderr)

    print("[CLADE] Stage 1: Fisher/FDR pre-filter...", file=sys.stderr)
    fisher = fisher_fdr_screen(df, candidates, pheno_col)
    sig = set(fisher.loc[fisher["FDR_significant"], "Candidate"])

    st_dummies = st_covariate_dummies(df)
    print("[CLADE] Stage 1: Firth regression for FDR-significant candidates...", file=sys.stderr)
    firth_results = {c: firth_association(df, c, pheno_col, st_dummies) for c in sig}

    print("[CLADE] Stage 2: lineage recurrence...", file=sys.stderr)
    recurrence = lineage_recurrence(df[candidates], df["ST"]) if "ST" in df.columns else None

    temporal = None
    if args.tree:
        print("[CLADE] Stage 4: temporal ordering (Fitch parsimony)...", file=sys.stderr)
        tree = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True)
        pheno_map = df[pheno_col].to_dict()
        cand_maps = {c: df[c].to_dict() for c in candidates}
        temporal = temporal_ordering(tree, pheno_map, cand_maps).set_index("Candidate")
        baseline = cohort_baseline_prevalence(pheno_map)
    else:
        print("[CLADE] Stage 4 skipped: no --tree supplied.", file=sys.stderr)
        baseline = None

    matched = None
    if args.distances:
        print("[CLADE] Stage 5: matched-neighbor comparison...", file=sys.stderr)
        dist = pd.read_csv(args.distances, sep="\t", index_col=0)
        dist.columns = dist.columns.astype(str)
        dist.index = dist.index.astype(str)
        pos = df[df[pheno_col] == 1].index.tolist()
        neg = df[df[pheno_col] == 0].index.tolist()
        common = [s for s in df.index if s in dist.index]
        pos = [s for s in pos if s in common]
        neg = [s for s in neg if s in common]
        nn = nearest_negative_neighbor(dist.loc[common, common], pos, neg)
        matched = matched_neighbor_test(df, pos, nn, candidates).set_index("Candidate")
    else:
        print("[CLADE] Stage 5 skipped: no --distances supplied.", file=sys.stderr)

    results = {}
    for c in candidates:
        firth = firth_results.get(c)
        stage1_tested = firth is not None
        stage1_sig = (firth is not None and firth["Status"] == "OK" and firth["Firth_p"] is not None and firth["Firth_p"] < 0.05)
        direction = direction_check(firth["Coef"]) if firth and firth["Status"] == "OK" else {"correct_direction": None}

        stage2_pct = recurrence.set_index("Candidate").loc[c, "Pct_dominant_ST"] if recurrence is not None and c in recurrence["Candidate"].values else None
        stage2_sts = recurrence.set_index("Candidate").loc[c, "Distinct_STs"] if recurrence is not None and c in recurrence["Candidate"].values else None

        pct_post = temporal.loc[c, "Pct_post_resistance"] if temporal is not None and c in temporal.index else None
        s5_p = matched.loc[c, "McNemar_p"] if matched is not None and c in matched.index else None
        s5_sig = (s5_p is not None and s5_p < 0.05)
        s5_concordant = None  # requires comparing direction of Stage 5 effect to Stage 3/4 - left for caller interpretation

        r = StageResults(
            stage1_tested=stage1_tested,
            stage1_significant=stage1_sig,
            stage2_distinct_sts=stage2_sts,
            stage2_pct_dominant_st=stage2_pct,
            stage3_correct_direction=direction["correct_direction"],
            stage4_pct_post_resistance=pct_post,
            stage4_baseline=baseline,
            stage5_significant=s5_sig,
            stage5_concordant=s5_concordant,
            cohort_frequency=float(df[c].sum()) / len(df),
        )
        results[c] = r

    table = build_evidence_table(results)
    markdown = evidence_table_to_markdown(table)

    with open(args.output, "w") as f:
        f.write(markdown)
    print(f"[CLADE] Wrote evidence table to {args.output}", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="clade", description="CLADE: Clonal-Lineage-Aware Detection of Epistasis")
    parser.add_argument("--version", action="version", version=f"clade {__version__}")
    subparsers = parser.add_subparsers(dest="command")

    p_validate = subparsers.add_parser("validate", help="Run the six-stage validation framework on a candidate list")
    p_validate.add_argument("--genotypes", nargs="+", required=True, help="One or more candidate-genotype CSVs")
    p_validate.add_argument("--phenotype", required=True, help="Phenotype TSV (sample -> 0/1)")
    p_validate.add_argument("--lineages", required=True, help="MLST TSV (sample, ST columns)")
    p_validate.add_argument("--candidates", nargs="+", required=True, help="Candidate column names to validate")
    p_validate.add_argument("--tree", help="Newick phylogeny (enables Stage 4)")
    p_validate.add_argument("--distances", help="Pairwise patristic distance TSV (enables Stage 5)")
    p_validate.add_argument("--config", help="Reserved for future use")
    p_validate.add_argument("--output", required=True, help="Output evidence-table Markdown path")
    p_validate.set_defaults(func=_cmd_validate)

    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
