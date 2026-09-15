"""
Generate the synthetic validation cohort in tests/fixtures/synthetic_cohort/.

**This is software-validation data, not biological evidence.** Nothing produced
here says anything about any organism. Its only purpose is to give the pipeline
inputs whose correct answer is known by construction, so a scientific
regression can be detected.

Four candidates are planted with deliberately different truths:

    true_compensatory  arises independently in three lineages, always on an
                       already-resistant background, enriched in resistant
                       genomes                      -> should reach CONVERGENT
    clonal_marker      confined to a single lineage, one ancestral origin,
                       inherited by every descendant -> should NOT be convergent
    wrong_direction    concentrated in susceptible genomes
                                                     -> explicit Stage 3 failure
    near_fixed         present in almost every genome
                                                     -> UNINFORMATIVE

Run:  python tests/fixtures/make_synthetic_cohort.py
Deterministic: fixed seed, no network, no external data.
"""
from __future__ import annotations

import json
import pathlib

import dendropy
import numpy as np
import pandas as pd

SEED = 20260913
OUT = pathlib.Path(__file__).parent / "synthetic_cohort"

# lineage -> (n_samples, n_resistant)
#
# Every lineage is deliberately *mixed*. Lineages that perfectly predict
# resistance make the lineage covariates separate the data completely, the
# logistic fit becomes non-identifiable, and the "structure-corrected" results
# below degenerate into absurd coefficients with unusable p-values. That is a
# real phenomenon (it is why this project needed Firth regression at all), but
# it is not what this fixture is for: here we want a cohort where a
# structure-corrected test can actually be computed.
LINEAGES = {"ST1": (16, 13), "ST2": (16, 12), "ST3": (16, 8), "ST4": (12, 2)}


def build() -> None:
    rng = np.random.default_rng(SEED)
    OUT.mkdir(exist_ok=True)

    samples, lineage_of, resistant = [], {}, {}
    for st, (n, n_res) in LINEAGES.items():
        for i in range(n):
            s = f"{st}_{i:02d}"
            samples.append(s)
            lineage_of[s] = st
            resistant[s] = 1 if i < n_res else 0

    # ---- phylogeny: each lineage is a clade, balanced within ----
    def clade(members: list[str]) -> str:
        if len(members) == 1:
            return f"{members[0]}:0.01"
        mid = len(members) // 2
        return f"({clade(members[:mid])},{clade(members[mid:])}):0.02"

    by_st = {st: [s for s in samples if lineage_of[s] == st] for st in LINEAGES}
    newick = (
        f"(({clade(by_st['ST1'])},{clade(by_st['ST2'])}):0.1,"
        f"({clade(by_st['ST3'])},{clade(by_st['ST4'])}):0.1);"
    )
    (OUT / "tree.nwk").write_text(newick + "\n")

    # ---- candidate genotypes ----
    geno = pd.DataFrame(0, index=samples, columns=[
        "true_compensatory", "clonal_marker", "wrong_direction", "near_fixed"
    ], dtype=int)
    geno.index.name = "Sample_ID"

    # true_compensatory: three independent origins, each inside a clade whose
    # ancestor is already resistant. A couple of susceptible carriers keep the
    # data from being perfectly separated (which would make any logistic fit
    # non-identifiable — the very problem Firth regression exists to handle).
    geno.loc[by_st["ST1"][:7], "true_compensatory"] = 1
    geno.loc[by_st["ST2"][:6], "true_compensatory"] = 1
    geno.loc[by_st["ST3"][:4], "true_compensatory"] = 1
    geno.loc[by_st["ST4"][-2:], "true_compensatory"] = 1  # susceptible carriers

    # clonal_marker: one ancestral gain, inherited by an entire lineage.
    geno.loc[by_st["ST1"], "clonal_marker"] = 1

    # wrong_direction: concentrated in the susceptible lineage.
    geno.loc[by_st["ST4"][2:], "wrong_direction"] = 1
    geno.loc[by_st["ST3"][8:], "wrong_direction"] = 1
    geno.loc[by_st["ST2"][-3:], "wrong_direction"] = 1
    geno.loc[by_st["ST1"][:2], "wrong_direction"] = 1

    # near_fixed: present almost everywhere.
    geno["near_fixed"] = 1
    geno.loc[rng.choice(samples, size=2, replace=False), "near_fixed"] = 0

    geno.to_csv(OUT / "genotypes.csv")

    pheno = pd.DataFrame({"sample": samples, "resistant": [resistant[s] for s in samples]})
    pheno.to_csv(OUT / "phenotype.tsv", sep="\t", index=False)

    mlst = pd.DataFrame({"sample": samples, "ST": [lineage_of[s] for s in samples]})
    mlst.to_csv(OUT / "lineages.tsv", sep="\t", index=False)

    # ---- real patristic distances from the tree we just wrote ----
    tree = dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)
    pdm = tree.phylogenetic_distance_matrix()
    label = {t.label: t for t in tree.taxon_namespace}
    dist = pd.DataFrame(
        [[pdm.distance(label[a], label[b]) if a != b else 0.0 for b in samples] for a in samples],
        index=samples,
        columns=samples,
    )
    dist.to_csv(OUT / "distances.tsv", sep="\t")

    # ---- structure-corrected Stage 1 results ----
    # Real logistic regressions with lineage indicator covariates, computed on
    # this synthetic data. They stand in for pyseer output so the preferred
    # Stage 1 path is exercised without requiring firthlogist.
    import statsmodels.api as sm

    st_dummies = pd.get_dummies(
        pd.Series([lineage_of[s] for s in samples], index=samples), prefix="ST", drop_first=True
    ).astype(float)
    y = pd.Series([resistant[s] for s in samples], index=samples).astype(int)

    rows = []
    for cand in geno.columns:
        X = sm.add_constant(pd.concat([st_dummies, geno[[cand]].astype(float)], axis=1))
        try:
            fit = sm.Logit(y, X).fit(disp=0, maxiter=200)
            coef, pval = float(fit.params[cand]), float(fit.pvalues[cand])
            if not np.isfinite(coef) or not np.isfinite(pval) or abs(coef) > 10:
                # Implausible magnitude means the fit did not identify: record
                # it as missing rather than shipping a meaningless number.
                coef, pval = float("nan"), float("nan")
        except Exception:  # noqa: BLE001 - any fit failure means 'no estimate'
            coef, pval = float("nan"), float("nan")
        rows.append({"Candidate": cand, "Coef": coef, "P": pval})
    pd.DataFrame(rows).to_csv(OUT / "stage1_results.tsv", sep="\t", index=False)

    # ---- recorded manual Stage 6 findings ----
    (OUT / "stage6.json").write_text(
        json.dumps(
            {"true_compensatory": True, "clonal_marker": None,
             "wrong_direction": None, "near_fixed": None},
            indent=2,
        )
        + "\n"
    )

    (OUT / "README.md").write_text(
        "# Synthetic validation cohort\n\n"
        "**Software validation data only — not biological evidence.**\n\n"
        f"Generated by `make_synthetic_cohort.py` (seed {SEED}); regenerate with\n"
        "`python tests/fixtures/make_synthetic_cohort.py`.\n\n"
        f"{len(samples)} samples across {len(LINEAGES)} lineages; "
        f"{int(sum(resistant.values()))} resistant.\n\n"
        "| Candidate | Planted truth | Expected disposition |\n"
        "|---|---|---|\n"
        "| `true_compensatory` | three independent origins, all post-resistance | convergent |\n"
        "| `clonal_marker` | one origin, inherited by one whole lineage | not convergent |\n"
        "| `wrong_direction` | concentrated in susceptible genomes | rejected |\n"
        "| `near_fixed` | present in almost every genome | uninformative |\n"
    )
    print(f"wrote {len(samples)} samples to {OUT}")


if __name__ == "__main__":
    build()
