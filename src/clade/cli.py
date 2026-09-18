"""
CLADE command-line interface.

    clade --version
    clade validate --genotypes G1.csv [G2.csv ...] --phenotype pheno.tsv \\
        --lineages mlst.tsv --candidates cand1 cand2 ... \\
        [--tree tree.nwk] [--distances distances.tsv] \\
        --output evidence_table.md

Stages 4 (temporal ordering) and 5 (matched-neighbour) are skipped, with a
clear message, if `--tree` / `--distances` are not supplied. A skipped stage is
recorded as *not run*, never as a negative result — a candidate cannot become
convergent because a stage was absent. Stage 6 (external corroboration) is
never run by the CLI: it is a manual literature/homology process by design (see
clade.validation.stage6_corroboration), and its findings are supplied through
`--stage6` if they have been recorded.

Exit codes:
    0  completed
    1  usage error (no subcommand)
    2  invalid or inconsistent input data
    3  a required optional dependency is missing
"""
from __future__ import annotations

import argparse
import json
import sys

from clade import __version__

EXIT_OK = 0
EXIT_USAGE = 1
EXIT_BAD_INPUT = 2
EXIT_MISSING_DEP = 3


def _log(msg: str) -> None:
    print(f"[CLADE] {msg}", file=sys.stderr)


def _load_external_stage1(path: str) -> dict[str, dict]:
    """Load structure-corrected Stage 1 results computed outside CLADE.

    CLADE's preferred Stage 1 method is kinship-corrected regression (pyseer),
    which CLADE does not reimplement. This reads those results back in so the
    preferred method can actually drive the framework, rather than the shipped
    Firth fallback being the only reachable path.

    Expected TSV columns: `Candidate`, `Coef` (effect size, signed), and one of
    `P` or `Q` (significance). Anything else is ignored.
    """
    import pandas as pd

    from clade.io.validation import CladeInputError, is_missing

    tab = pd.read_csv(path, sep="\t")
    if "Candidate" not in tab.columns:
        raise CladeInputError(f"{path} has no 'Candidate' column; found {list(tab.columns)[:8]}")
    sig_col = next((c for c in ("Q", "q", "P", "p") if c in tab.columns), None)
    if sig_col is None:
        raise CladeInputError(f"{path} needs a 'P' or 'Q' column; found {list(tab.columns)[:8]}")
    if "Coef" not in tab.columns:
        raise CladeInputError(
            f"{path} needs a signed 'Coef' column — Stage 3 reads its direction from it."
        )

    out = {}
    for row in tab.itertuples(index=False):
        coef = row.Coef
        pval = getattr(row, sig_col)
        ok = not is_missing(coef) and not is_missing(pval)
        out[str(row.Candidate)] = {
            "Candidate": str(row.Candidate),
            "Coef": float(coef) if ok else None,
            "Firth_p": float(pval) if ok else None,
            "Status": "OK" if ok else "FAILED: missing estimate in external results",
        }
    return out


def _cmd_validate(args: argparse.Namespace) -> int:
    import dendropy
    import pandas as pd

    from clade.classification.disposition import StageResults
    from clade.io.genotypes import load_joined_dataset, st_covariate_dummies
    from clade.io.validation import (
        CladeInputError,
        is_missing,
        validate_distance_matrix,
        validate_genotype_frame,
        validate_phenotype,
    )
    from clade.reporting.evidence_table import build_evidence_table, evidence_table_to_markdown
    from clade.validation.stage1_association import firth_association, fisher_fdr_screen
    from clade.validation.stage2_recurrence import lineage_recurrence
    from clade.validation.stage3_direction import direction_check
    from clade.validation.stage4_temporal_order import (
        cohort_baseline_prevalence,
        temporal_ordering,
    )
    from clade.validation.stage5_matched_neighbors import (
        matched_neighbor_test,
        nearest_negative_neighbor,
    )

    pheno_col = "phenotype"
    candidates = list(dict.fromkeys(args.candidates))
    if len(candidates) != len(args.candidates):
        _log("note: duplicate candidate names collapsed.")

    df = load_joined_dataset(args.phenotype, args.genotypes, args.lineages)
    if df.empty:
        raise CladeInputError(
            "the phenotype, genotype and lineage files share no sample IDs. "
            "Check that all three use the same identifier convention."
        )

    validate_genotype_frame(df, candidates)
    validate_phenotype(df[pheno_col])
    _log(f"Loaded {len(df)} samples, {len(candidates)} candidates.")

    # ---- Stage 1 -------------------------------------------------------
    _log("Stage 1: Fisher/FDR prevalence screen...")
    screen = fisher_fdr_screen(
        df,
        candidates,
        pheno_col,
        min_carriers=args.min_carriers,
        max_freq=args.max_freq,
        alpha=args.alpha,
    ).set_index("Candidate")

    for status, n in screen["Screen_status"].value_counts().items():
        _log(f"  screen: {n} candidate(s) {status}")

    # The Fisher/FDR screen is NOT structure-corrected and can never, on its
    # own, satisfy Stage 1 — that is the entire point of the stage. It only
    # selects which candidates are worth a structure-corrected test.
    sig = [c for c in candidates if bool(screen.loc[c, "FDR_significant"])]

    external: dict[str, dict] = {}
    if args.stage1_results:
        external = _load_external_stage1(args.stage1_results)
        _log(
            f"Stage 1: loaded {len(external)} externally-computed, structure-corrected "
            f"result(s) from {args.stage1_results}"
        )

    firth_results: dict[str, dict] = dict(external)
    to_fit = [c for c in sig if c not in firth_results]
    if to_fit:
        _log(f"Stage 1: Firth regression for {len(to_fit)} screen-significant candidate(s)...")
        st_dummies = st_covariate_dummies(df)
        try:
            for c in to_fit:
                firth_results[c] = firth_association(df, c, pheno_col, st_dummies)
        except ImportError:
            _log(
                "WARNING: Stage 1's Firth fallback needs the optional 'firthlogist' package, "
                "which is not installed (it requires Python <3.11 and scikit-learn <1.6).\n"
                "         Install with: pip install 'clade[firth]'\n"
                "         Or supply kinship-corrected results from pyseer via --stage1-results.\n"
                "         Stage 1 is recorded as INDETERMINATE for every candidate: the "
                "uncorrected Fisher/FDR screen is a pre-filter, not a substitute for a "
                "structure-corrected test, and CLADE will not report it as one."
            )
            for c in to_fit:
                firth_results[c] = {
                    "Candidate": c,
                    "Coef": None,
                    "Firth_p": None,
                    "Status": "FAILED: firthlogist not installed",
                }
    elif not sig:
        _log("Stage 1: no candidate passed the screen into structure-corrected testing.")

    # ---- Stage 2 -------------------------------------------------------
    recurrence = None
    if "ST" in df.columns:
        _log("Stage 2: lineage recurrence...")
        recurrence = lineage_recurrence(df[candidates], df["ST"]).set_index("Candidate")
        warning = recurrence.attrs.get("lineage_warning")
        if warning:
            _log(f"  warning: {warning}")
    else:
        _log("Stage 2 skipped: no ST column available.")

    # ---- Stage 4 -------------------------------------------------------
    temporal = None
    baseline = None
    if args.tree:
        _log("Stage 4: temporal ordering (Fitch parsimony)...")
        tree = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True)
        pheno_map = df[pheno_col].to_dict()
        cand_maps = {c: df[c].to_dict() for c in candidates}
        temporal = temporal_ordering(
            tree,
            pheno_map,
            cand_maps,
            ambiguity_resolution=args.ambiguity_resolution,
            min_overlap_fraction=args.min_tree_overlap,
        ).set_index("Candidate")
        baseline = cohort_baseline_prevalence(pheno_map)
        _log(f"  cohort resistance baseline: {baseline:.3f}")
    else:
        _log("Stage 4 skipped: no --tree supplied (recorded as not run, not as a failure).")

    # ---- Stage 5 -------------------------------------------------------
    matched = None
    if args.distances:
        _log("Stage 5: matched-neighbour comparison...")
        dist = pd.read_csv(args.distances, sep="\t", index_col=0)
        dist.index = dist.index.astype(str)
        dist.columns = dist.columns.astype(str)
        dist = validate_distance_matrix(dist)

        in_matrix = [s for s in df.index if s in set(dist.index)]
        if not in_matrix:
            raise CladeInputError(
                "no sample in the dataset appears in the distance matrix; "
                "check the identifier convention."
            )
        if len(in_matrix) < len(df):
            _log(f"  note: {len(df) - len(in_matrix)} sample(s) absent from the distance matrix.")

        sub = df.loc[in_matrix]
        pos = sub.index[sub[pheno_col] == 1].tolist()
        neg = sub.index[sub[pheno_col] == 0].tolist()
        if not pos or not neg:
            _log("  Stage 5 skipped: need both resistant and susceptible samples to form pairs.")
        else:
            nn = nearest_negative_neighbor(
                dist.loc[in_matrix, in_matrix], pos, neg, unique=args.unique_neighbors
            )
            matched = matched_neighbor_test(df, pos, nn, candidates).set_index("Candidate")
            reuse = int(matched["Max_neighbor_reuse"].iloc[0])
            if reuse > 1:
                _log(
                    f"  warning: one susceptible genome anchors up to {reuse} pairs "
                    "(pseudoreplication; McNemar assumes independent pairs). "
                    "Re-run with --unique-neighbors as a sensitivity check."
                )
    else:
        _log("Stage 5 skipped: no --distances supplied (recorded as not run, not as a failure).")

    # ---- Stage 6 (manual, supplied) ------------------------------------
    stage6: dict[str, bool | None] = {}
    if args.stage6:
        with open(args.stage6) as fh:
            stage6 = json.load(fh)
        _log(f"Stage 6: loaded {len(stage6)} recorded corroboration result(s).")

    # ---- assemble ------------------------------------------------------
    results = {}
    for c in candidates:
        screen_status = str(screen.loc[c, "Screen_status"])
        firth = firth_results.get(c)

        # Stage 1 bookkeeping, kept deliberately explicit. The three states are
        # not interchangeable and collapsing any pair of them was the original
        # defect:
        #
        #   excluded by the prevalence screen -> never tested   (NOT_RETESTED)
        #   screened, no raw association      -> tested, False  (REJECTED)
        #   screened, association, no valid
        #     structure-corrected estimate    -> tested, None   (INSUFFICIENT)
        #   structure-corrected estimate      -> tested, True/False
        #
        # The Fisher/FDR screen is uncorrected, so it may only ever produce a
        # negative Stage 1 result, never a positive one: structure correction
        # can weaken an association but not create one.
        if screen_status != "tested":
            stage1_tested = False
            stage1_sig = None
        elif firth is None:
            # Screened and found to have no raw association at all.
            stage1_tested = True
            stage1_sig = False
        elif firth["Firth_p"] is not None:
            # A p-value is present exactly when the design is identifiable
            # enough to report one -- `firth_association` sets Firth_p to
            # None only for complete separation (coefficient undefined).
            # Quasi-separation (a handful of cases in the minority cell)
            # still yields a real, if less stable, p-value and is not
            # collapsed into "no estimate": it is reported, flagged unstable
            # in Status, and left to the reader with that caveat -- exactly
            # the distinction HTZ92_2925 M21L needed (min cell 2, otherwise
            # a clean Stage 1 pass) but the aminotransferase candidate did
            # not get (min cell 0, genuinely no estimate).
            stage1_tested = True
            stage1_sig = bool(firth["Firth_p"] < args.alpha)
        else:
            # A structure-corrected test was attempted but yielded no estimate
            # (complete separation, or the fit itself failed to converge).
            stage1_tested = True
            stage1_sig = None

        direction = (
            direction_check(firth["Coef"])
            if firth is not None and firth["Coef"] is not None
            else {"correct_direction": None}
        )

        s2_pct = s2_sts = None
        if recurrence is not None and c in recurrence.index:
            s2_pct = recurrence.loc[c, "Pct_dominant_ST"]
            s2_sts = recurrence.loc[c, "Distinct_STs"]

        pct_post = total_gains = None
        if temporal is not None and c in temporal.index:
            pct_post = float(temporal.loc[c, "Pct_post_resistance"])
            total_gains = int(temporal.loc[c, "Total_gains"])

        s5_sig = s5_dir = None
        if matched is not None and c in matched.index:
            p = matched.loc[c, "McNemar_p"]
            if not is_missing(p):
                s5_sig = bool(p < args.alpha)
                s5_dir = matched.loc[c, "Enriched_in_resistant"]
                if s5_dir is not None:
                    s5_dir = bool(s5_dir)

        col = df[c]
        n_obs = int(col.notna().sum())
        results[c] = StageResults(
            stage1_tested=stage1_tested,
            stage1_significant=stage1_sig,
            stage2_distinct_sts=None if s2_sts is None else int(s2_sts),
            stage2_pct_dominant_st=None if s2_pct is None else float(s2_pct),
            stage2_independent_origins=total_gains,
            stage3_correct_direction=direction["correct_direction"],
            stage4_pct_post_resistance=pct_post,
            stage4_baseline=baseline,
            stage4_total_gains=total_gains,
            stage5_significant=s5_sig,
            stage5_enriched_in_resistant=s5_dir,
            stage6_mechanism_found=stage6.get(c),
            cohort_frequency=(float(col.sum()) / n_obs) if n_obs else None,
        )

    table = build_evidence_table(results, temporal_margin=args.temporal_margin)
    with open(args.output, "w") as f:
        f.write(evidence_table_to_markdown(table))
    _log(f"Wrote evidence table to {args.output}")

    if args.output_tsv:
        table.to_csv(args.output_tsv, sep="\t", index=False)
        _log(f"Wrote machine-readable table to {args.output_tsv}")

    for tier, n in table["Disposition"].value_counts().items():
        _log(f"  {tier}: {n}")
    return EXIT_OK


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="clade",
        description="CLADE: Clonal-Lineage-Aware Detection of Epistasis — "
        "six-stage validation of candidate compensatory mutations.",
    )
    parser.add_argument("--version", action="version", version=f"clade {__version__}")
    subparsers = parser.add_subparsers(dest="command")

    p = subparsers.add_parser(
        "validate",
        help="Run the six-stage validation framework on a candidate list",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    req = p.add_argument_group("required inputs")
    req.add_argument("--genotypes", nargs="+", required=True, help="One or more candidate-genotype CSVs")
    req.add_argument("--phenotype", required=True, help="Phenotype TSV (sample -> 0/1)")
    req.add_argument("--lineages", required=True, help="MLST TSV with sample and ST columns")
    req.add_argument("--candidates", nargs="+", required=True, help="Candidate column names to validate")
    req.add_argument("--output", required=True, help="Output evidence-table Markdown path")

    opt = p.add_argument_group("optional stage inputs")
    opt.add_argument("--tree", help="Newick phylogeny (enables Stage 4)")
    opt.add_argument("--distances", help="Pairwise patristic distance TSV (enables Stage 5)")
    opt.add_argument(
        "--stage1-results",
        help="TSV of externally-computed structure-corrected Stage 1 results "
        "(columns: Candidate, Coef, and P or Q) — e.g. pyseer output. This is "
        "CLADE's preferred Stage 1 path; the built-in Firth fit is the fallback.",
    )
    opt.add_argument("--stage6", help="JSON file: {candidate: true|false|null} of manual Stage 6 findings")
    opt.add_argument("--output-tsv", help="Also write the full table as TSV")

    thr = p.add_argument_group("thresholds")
    thr.add_argument("--alpha", type=float, default=0.05, help="Significance level")
    thr.add_argument("--min-carriers", type=int, default=10, help="Stage 1 screen: minimum carriers")
    thr.add_argument("--max-freq", type=float, default=0.80, help="Stage 1 screen: maximum carrier frequency")
    thr.add_argument(
        "--temporal-margin",
        type=float,
        default=0.10,
        help="Stage 4: how far Pct_post_resistance must exceed the cohort baseline",
    )
    thr.add_argument(
        "--ambiguity-resolution",
        type=int,
        choices=(0, 1),
        default=0,
        help="Fitch tie-breaking: 0 dates gains late (DELTRAN-like), 1 early (ACCTRAN-like)",
    )
    thr.add_argument(
        "--min-tree-overlap",
        type=float,
        default=0.5,
        help="Stage 4: minimum fraction of samples that must be present as tree leaves",
    )
    thr.add_argument(
        "--unique-neighbors",
        action="store_true",
        help="Stage 5: greedy one-to-one matching instead of allowing neighbour reuse "
        "(sensitivity check against pseudoreplication)",
    )
    p.set_defaults(func=_cmd_validate)

    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return EXIT_USAGE

    from clade.io.validation import CladeInputError

    try:
        return args.func(args)
    except CladeInputError as e:
        print(f"[CLADE] input error: {e}", file=sys.stderr)
        return EXIT_BAD_INPUT
    except FileNotFoundError as e:
        print(f"[CLADE] file not found: {e.filename}", file=sys.stderr)
        return EXIT_BAD_INPUT


if __name__ == "__main__":
    sys.exit(main())
