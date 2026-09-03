#!/usr/bin/env python3
"""
CAPSTONE · SCRIPT 3 of 3 — Quality-Assurance Audit: Compare Against the
Published Risk-Score Coefficients
==============================================================================
NEW SCRIPT. This is "Phase 6: Quality-assurance audit" of this capstone's 8
phases, made concrete: the source publication (Xu, Zheng & Xu 2026, Sci Rep,
https://doi.org/10.1038/s41598-026-68769-z - see this project's README for
the full citation and how it was located) reports the exact 9 gene
coefficients its own LASSO Cox model converged to (its Fig. 2A/2C). Re-running
PRJNA482620_TCGA_BRCA_downstream.R's Phase E (LASSO Cox regression) against a
freshly-pulled TCGA-BRCA cohort (Script 1) will NOT reproduce those exact
numbers - GDC's live database has drifted since the paper's own download date
(see Script 1's header), LASSO's cross-validation fold assignment is
stochastic unless seeded identically, and this project deliberately does not
replicate the paper's exact sample-filtering criteria end-to-end. The
question this script actually answers is narrower and more useful for a
capstone QA phase: do the SIGNS and RELATIVE MAGNITUDES of the re-derived
coefficients agree with the published model, or does a re-run produce a
qualitatively different signature - which would be a real finding worth
investigating, not something to paper over by re-running until the numbers
look right.

TEACHING NOTE — why sign agreement matters more than exact coefficient
match here:
  A LASSO Cox coefficient's SIGN says whether higher expression of that gene
  is associated with worse (positive) or better (negative) survival in this
  cohort - that is the actual scientific claim the published model makes
  about each gene. Its magnitude is comparatively fragile (rescales with
  cohort size, exact TPM normalization, and the regularization path's own
  cross-validation fold assignment). A coefficient disagreeing in SIGN from
  the published model is a real discrepancy worth investigating (data
  version drift, a sample-inclusion criterion difference, or a genuine
  script bug); a difference in magnitude alone, with the same sign, is
  expected variation, not evidence of a problem.

Input:  --derived-tsv (required): a TSV with columns 'gene' and
        'coefficient', produced by exporting PRJNA482620_TCGA_BRCA_downstream.R's
        cv_fit object, e.g.:
            coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))
            write.table(data.frame(gene = rownames(coefs), coefficient = coefs[, 1]),
                        "derived_coefficients.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
        (PRJNA482620_TCGA_BRCA_downstream.R does not write this file itself -
        exporting it is a real, disclosed gap in that script, not something
        silently patched here; add the snippet above before running Script 3
        of PRJNA482620_TCGA_BRCA_downstream.R for a real comparison.)
Output: prints a per-gene sign/magnitude comparison; writes
        coefficient_qa_summary.tsv
"""

import argparse
import sys
from pathlib import Path

# Published coefficients (Xu, Zheng & Xu 2026, Fig. 2C) - the paper's own
# 9-gene LASSO Cox signature, transcribed directly from the article text.
PUBLISHED_COEFFICIENTS = {
    "CCDC74B": -0.095,
    "ELOVL2": -0.070,
    "FAM234B": 0.248,
    "FUT3": 0.028,
    "KCNJ11": -0.041,
    "NXPH4": 0.054,
    "RIBC1": -0.084,
    "UPK2": 0.101,
    "WNK4": -0.016,
}


def read_derived_tsv(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}
    for col in ("gene", "coefficient"):
        if col not in idx:
            raise ValueError(f"expected column '{col}' not found in {path} header: {header}")
    out = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        gene = cols[idx["gene"]]
        try:
            coef = float(cols[idx["coefficient"]])
        except ValueError:
            continue  # e.g. glmnet's "(Intercept)" row, or a zeroed-out gene
        if coef != 0.0:
            out[gene] = coef
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Compare re-derived LASSO Cox coefficients against the published 9-gene signature.")
    parser.add_argument("--derived-tsv", required=True,
                         help="derived_coefficients.tsv (gene, coefficient columns) from a real re-run")
    parser.add_argument("--out", default="coefficient_qa_summary.tsv")
    args = parser.parse_args()

    print("==================================================================")
    print(" CAPSTONE / SCRIPT 3 — Coefficient QA audit vs. published signature")
    print("==================================================================")

    if not Path(args.derived_tsv).exists():
        print(f"[FAIL] {args.derived_tsv} not found — export it from a real Script 3 "
              f"(PRJNA482620_TCGA_BRCA_downstream.R) run first; see this script's header for the exact snippet.")
        sys.exit(1)

    derived = read_derived_tsv(args.derived_tsv)
    all_genes = sorted(set(PUBLISHED_COEFFICIENTS) | set(derived))

    rows = []
    n_sign_disagree = 0
    n_published_missing = 0
    n_derived_extra = 0
    for gene in all_genes:
        pub = PUBLISHED_COEFFICIENTS.get(gene)
        der = derived.get(gene)
        if pub is None:
            n_derived_extra += 1
            print(f"  {gene:10s} published=NOT_IN_PUBLISHED_SIGNATURE  derived={der:+.3f}   NEW_GENE")
            rows.append((gene, "NA", f"{der:.3f}", "NEW_GENE"))
            continue
        if der is None:
            n_published_missing += 1
            print(f"  {gene:10s} published={pub:+.3f}  derived=DROPPED_BY_LASSO   MISSING_FROM_RERUN")
            rows.append((gene, f"{pub:.3f}", "NA", "MISSING_FROM_RERUN"))
            continue
        same_sign = (pub > 0) == (der > 0)
        flag = "ok" if same_sign else "SIGN_DISAGREEMENT"
        if not same_sign:
            n_sign_disagree += 1
        print(f"  {gene:10s} published={pub:+.3f}  derived={der:+.3f}   {flag}")
        rows.append((gene, f"{pub:.3f}", f"{der:.3f}", flag))

    print()
    if n_sign_disagree > 0:
        print(f"[WARN] {n_sign_disagree} gene(s) disagree in sign with the published model — "
              f"investigate before trusting this re-run's risk stratification, don't average it away.")
    else:
        print("[PASS] no gene disagrees in sign with the published model.")
    if n_published_missing > 0:
        print(f"[WARN] {n_published_missing} published gene(s) were dropped entirely by this re-run's LASSO "
              f"(coefficient shrunk to exactly 0) — expected occasionally with cohort drift, worth noting "
              f"in the capstone's Interpretation section either way.")
    if n_derived_extra > 0:
        print(f"[WARN] {n_derived_extra} gene(s) in this re-run were NOT in the published 9-gene signature — "
              f"a re-run legitimately can converge on a different gene set; this is not itself an error.")

    out_path = Path(args.out)
    with open(out_path, "w") as fh:
        fh.write("gene\tpublished_coefficient\tderived_coefficient\tflag\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    print(f"\nWrote {out_path}")
    print("==================================================================")


if __name__ == "__main__":
    main()
