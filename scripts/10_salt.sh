#!/bin/bash
# Step 10. Salt response (Methods and Results: "Salt response"; Figure S1, Table S5). Needs steps 1 and 3.
#  10.1 b_salt_B.R             (1) salt response of B genes vs donors and host genes; (2) carrier x salt on host genes
#  10.2 b_salt_B.R PERM_ONLY   1,000 permutations of (2): carrier labels shuffled among genotypes within lineages
#                              (100 seeds x 10; each seed can run as a separate job) -> B_salt/perm/
#  10.3 salt_perm_summary.R    observed vs permuted -> B_salt/perm_summary.tsv
#  10.4 b_salt_Btype.R         (3) salt response of B genes by B type, with host and donor controls
source "$(dirname "$0")/../config.sh"; set -euo pipefail
mkdir -p B_salt/perm
echo "== 10.1"; Rscript $S/b_salt_B.R
echo "== 10.2"
seq 1 100 | xargs -P $PAR -I{} sh -c 'PERM_ONLY=1 PERM_SEED={} NPERM=10 Rscript $S/b_salt_B.R > B_salt/perm/log_{}.txt 2>&1 || echo "seed {} failed"'
echo "== 10.3"; Rscript $S/salt_perm_summary.R
echo "== 10.4"; NPERM=0 Rscript $S/b_salt_Btype.R
echo "=== step 10 done"
