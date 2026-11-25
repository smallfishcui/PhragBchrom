#!/bin/bash
# Step 4. GO enrichment of B genes, counted per gene and per B gene family; leave-one-family-out
#         (Methods: "Functional enrichment"; Table S3, Figure 3b). Needs step 3.
#  4.1 b_go_main.R   -> Bgenes/GO_main_gene.tsv, Bgenes/GO_main_family.tsv
source "$(dirname "$0")/../config.sh"; set -euo pipefail
mkdir -p Bgenes
echo "== 4.1"; Rscript $S/b_go_main.R
echo "=== step 4 done"
