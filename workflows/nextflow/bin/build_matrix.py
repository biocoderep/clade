#!/usr/bin/env python3
"""
Merge per-sample outputs into CLADE's three required inputs: a candidate
genotype matrix, a phenotype vector, and a lineage file.

Genuinely new script. Two real design points worth being explicit about:

1. Candidate genotypes come from `core.vcf` by exact chromosome:position
   match, NOT from Abricate/AMRFinder/MLST tables — those give gene-level
   presence/absence, useful for the phenotype, but not the variant-level
   calls CLADE's candidates need (e.g. "SecA_M21L" is a specific amino-acid
   substitution, not just "is SecA present"). This mirrors exactly how the
   real case study's `extract_remaining_candidates.py` worked.

2. The phenotype is derived, not hand-maintained, reusing
   `clade.provenance.phenotype.derive_phenotype_from_amr_matrix` — the same
   function that encodes the fix for a real phenotype-identity error caught
   during the case study (see docs/reproducibility/phenotype_correction.md).
   This script's only job for the phenotype is to aggregate per-sample
   AMRFinder tables into the Rtab-style matrix that function expects.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "src"))
from clade.provenance.phenotype import derive_phenotype_from_amr_matrix


def parse_candidates(path: str) -> list[tuple[str, str, int]]:
    """Read a `name,position` CSV (position = CHROM:POS) into (name, chrom, pos) tuples."""
    out = []
    with open(path) as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].strip().lower() == "name":
                continue
            name, position = row[0].strip(), row[1].strip()
            chrom, pos = position.split(":")
            out.append((name, chrom, int(pos)))
    return out


#: VCF genotype strings meaning "reference allele called here".
REF_CALLS = frozenset({"0", "0/0", "0|0"})
#: VCF genotype strings meaning "no call" — the position could not be
#: genotyped for this sample. These are NOT reference calls.
NO_CALLS = frozenset({".", "./.", ".|."})


def extract_genotypes(vcf_path: str, candidates: list[tuple[str, str, int]]) -> pd.DataFrame:
    """Extract 0/1 genotype calls at each candidate position from a (multi-sample) VCF.

    Missing data is preserved as missing, not coerced to absence. Two cases
    matter, and both used to silently become 0:

    * a **no-call** (`.`, `./.`) means the position could not be genotyped in
      that sample — usually low coverage. Recording it as 0 asserts the sample
      does not carry the mutation, which the data does not support, and
      systematically inflates the non-carrier group.
    * a candidate position **absent from the VCF entirely** (filtered, or
      outside the core genome) yields no information about any sample. The
      whole column is missing, not all-zero.

    CLADE's stages treat a missing genotype as missing throughout: Stage 4
    reconstructs it as an ambiguous ancestral state, Stage 5 drops the pair,
    and Stage 1 drops the sample from that candidate's fit.
    """
    with open(vcf_path) as f:
        samples: list[str] = []
        rows: dict[str, dict[str, float]] = {name: {} for name, _, _ in candidates}
        wanted = {(chrom, pos): name for name, chrom, pos in candidates}
        seen_positions: set[str] = set()

        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                # Start every candidate as entirely missing. A position never
                # seen in the VCF stays missing rather than becoming all-zero.
                for name in rows:
                    rows[name] = {s: float("nan") for s in samples}
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos = fields[0], int(fields[1])
            key = (chrom, pos)
            if key not in wanted:
                continue
            name = wanted[key]
            seen_positions.add(name)
            genotype_fields = fields[9:]
            for sample, gt_field in zip(samples, genotype_fields):
                gt = gt_field.split(":")[0]
                if gt in NO_CALLS:
                    rows[name][sample] = float("nan")
                elif gt in REF_CALLS:
                    rows[name][sample] = 0.0
                else:
                    rows[name][sample] = 1.0

    absent = sorted(set(rows) - seen_positions)
    if absent:
        print(
            f"WARNING: {len(absent)} candidate position(s) not present in {vcf_path}: "
            f"{absent[:5]}{' ...' if len(absent) > 5 else ''}. "
            "Their genotypes are recorded as missing, not as absent.",
            file=sys.stderr,
        )

    df = pd.DataFrame(rows)
    df.index.name = "Sample_ID"
    return df.reset_index()


def build_amr_rtab(profiling_dir: str) -> pd.DataFrame:
    """Aggregate per-sample *_amrfinder.tsv files into an Rtab-style matrix
    (rows = gene symbol, columns = sample, first column header = 'Gene'),
    matching the format derive_phenotype_from_amr_matrix expects."""
    gene_presence: dict[str, dict[str, int]] = {}
    samples: list[str] = []

    for path in sorted(glob.glob(os.path.join(profiling_dir, "*_amrfinder.tsv"))):
        sample = os.path.basename(path).replace("_amrfinder.tsv", "")
        samples.append(sample)
        try:
            table = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            table = pd.DataFrame(columns=["Gene symbol"])
        genes_present = set(table.get("Gene symbol", pd.Series(dtype=str)).dropna())
        for gene in genes_present:
            gene_presence.setdefault(gene, {})[sample] = 1

    all_genes = sorted(gene_presence)
    rows = []
    for gene in all_genes:
        row = {"Gene": gene}
        for sample in samples:
            row[sample] = gene_presence[gene].get(sample, 0)
        rows.append(row)
    return pd.DataFrame(rows, columns=["Gene"] + samples)


def build_lineages(profiling_dir: str) -> pd.DataFrame:
    """Aggregate per-sample *_mlst.tsv files into a (sample, ST) table."""
    rows = []
    for path in sorted(glob.glob(os.path.join(profiling_dir, "*_mlst.tsv"))):
        sample = os.path.basename(path).replace("_mlst.tsv", "")
        with open(path) as f:
            line = f.readline().strip()
        st = None
        if line:
            parts = line.split("\t")
            if len(parts) >= 3:
                st = parts[2]
        rows.append({"sample": sample, "ST": st if st else "unknown"})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="core.vcf from SNIPPY_CORE")
    parser.add_argument("--candidates-file", required=True, help="CSV: name,position (CHROM:POS)")
    parser.add_argument("--profiling-dir", required=True, help="Directory with *_amrfinder.tsv, *_mlst.tsv")
    parser.add_argument("--resistance-gene", required=True, help="AMRFinderPlus gene symbol defining the phenotype")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    candidates = parse_candidates(args.candidates_file)
    genotypes = extract_genotypes(args.vcf, candidates)
    genotypes.to_csv(os.path.join(args.outdir, "genotypes.csv"), index=False)
    print(f"Wrote genotype matrix: {len(genotypes)} samples x {len(candidates)} candidates")

    amr_rtab_path = os.path.join(args.outdir, "amr_matrix.Rtab")
    build_amr_rtab(args.profiling_dir).to_csv(amr_rtab_path, sep="\t", index=False)
    phenotype = derive_phenotype_from_amr_matrix(amr_rtab_path, args.resistance_gene)
    phenotype.to_frame(name="phenotype").to_csv(os.path.join(args.outdir, "pheno.tsv"), sep="\t")
    print(f"Wrote phenotype vector: {int(phenotype.sum())}/{len(phenotype)} positive for {args.resistance_gene}")

    lineages = build_lineages(args.profiling_dir)
    lineages.to_csv(os.path.join(args.outdir, "lineages.tsv"), sep="\t", index=False)
    print(f"Wrote lineage file: {len(lineages)} samples")


if __name__ == "__main__":
    main()
