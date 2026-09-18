#!/usr/bin/env python3
"""Generate a genome-wide candidates.csv directly from core.vcf, replacing
the requirement for a user-supplied candidate list.

Why this exists, and what it deliberately does NOT do: this project's own
original discovery step (Mutual Information + Direct Coupling Analysis on a
synthetic sequence encoding) was found invalid -- no population-structure
correction at all, so any variant fixed early in a successful clonal lineage
scored as a top hit regardless of function (see the manuscript's §5.5 and
§7.6). This script does not repeat that mistake by trying to be a smarter
pre-ranking step. It does the opposite: it makes NO ranking claim whatsoever.
Every biallelic SNP in the core alignment becomes one candidate, unranked,
and the entire burden of separating signal from clonal noise is left to
CLADE's six stages -- Fisher/FDR prevalence screening, then structure-
corrected association, then the four checks that generate the ranking a
discovery step should never have made on its own.

Practical consequence: this is usually a LARGE candidate list (one row per
core-genome SNP -- tens of thousands on a real cohort, not the curated
dozens this framework's case study worked from). That is by design, not a
bug to fix by pre-filtering here. `fisher_fdr_screen` inside CLADE Stage 1
already does the prevalence-based pre-filtering (min_carriers / max_freq),
which is where structure-blind filtering belongs, not here, where it would
run before any structure correction exists to check it against.

What IS filtered out here, and why it is safe to do so before any
structure-aware step: monomorphic sites (no information in any test),
multi-allelic sites (CLADE's genotype loader requires strictly binary
columns and a triallelic SNP has no single meaningful 0/1 encoding), and
sites with excessive missingness (--max-missing-fraction), because a site
where most samples have no call cannot inform any test regardless of
population structure -- this is a data-completeness filter, not a
signal-based one.

Usage:
    discover_candidates.py --vcf core.vcf --output candidates.csv \
        --max-missing-fraction 0.10
"""
from __future__ import annotations

import argparse
import sys

REF_CALLS = frozenset({"0", "0/0", "0|0"})
ALT_CALLS = frozenset({"1", "1/1", "1|1"})
NO_CALLS = frozenset({".", "./.", ".|."})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="core.vcf from SNIPPY_CORE")
    ap.add_argument("--output", required=True, help="candidates.csv: name,position")
    ap.add_argument("--max-missing-fraction", type=float, default=0.10,
                    help="drop sites with more than this fraction of no-calls "
                         "across the cohort (default 0.10)")
    ap.add_argument("--max-sites", type=int, default=None,
                    help="hard cap on candidates emitted, for a quick smoke test; "
                         "sites are kept in VCF order, not by any score, since this "
                         "script makes no ranking claim -- omit for a real run")
    args = ap.parse_args()

    n_seen = n_kept = n_monomorphic = n_multiallelic = n_too_missing = 0
    rows = []

    with open(args.vcf) as fh:
        n_samples = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                n_samples = len(line.rstrip("\n").split("\t")) - 9
                continue
            if n_samples is None:
                sys.exit(f"{args.vcf}: no #CHROM header found before the first data line")

            n_seen += 1
            parts = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = parts[0:5]

            if "," in alt or len(alt) != 1 or len(ref) != 1:
                n_multiallelic += 1
                continue

            calls = parts[9:]
            n_ref = n_alt = n_no = 0
            for c in calls:
                gt = c.split(":", 1)[0].strip()
                if gt in NO_CALLS:
                    n_no += 1
                elif gt in REF_CALLS:
                    n_ref += 1
                elif gt in ALT_CALLS:
                    n_alt += 1
                else:
                    n_no += 1  # an unrecognised call is treated as missing, not guessed at

            if n_ref == 0 or n_alt == 0:
                n_monomorphic += 1
                continue
            if n_samples and (n_no / n_samples) > args.max_missing_fraction:
                n_too_missing += 1
                continue

            rows.append((f"{chrom}_{pos}", f"{chrom}:{pos}"))
            n_kept += 1
            if args.max_sites and n_kept >= args.max_sites:
                break

    if not rows:
        sys.exit(
            f"No candidate sites survived filtering out of {n_seen} core-genome variants "
            f"({n_monomorphic} monomorphic, {n_multiallelic} multiallelic, "
            f"{n_too_missing} too much missing data). Check --vcf is really a variant-site "
            f"VCF (snippy-core's core.vcf, not core.full.aln) and that "
            f"--max-missing-fraction is not too strict for this cohort."
        )

    with open(args.output, "w") as fh:
        fh.write("name,position\n")
        for name, pos in rows:
            fh.write(f"{name},{pos}\n")

    print(f"core.vcf: {n_seen} total variant sites, {n_samples} samples")
    print(f"  dropped monomorphic  : {n_monomorphic}")
    print(f"  dropped multiallelic : {n_multiallelic}")
    print(f"  dropped >{args.max_missing_fraction:.0%} missing : {n_too_missing}")
    print(f"  candidates emitted   : {n_kept}")
    print(f"Wrote {args.output}. This list carries no ranking -- CLADE's own six stages "
          f"are what separate signal from clonal noise, not this step.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
