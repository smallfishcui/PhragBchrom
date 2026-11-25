#!/bin/bash
# Step 7. B types and their divergence relative to the host (Methods: "B types and divergence"; Table 4, Figure 7,
#         Table S6). Needs step 6.
#  7.1 b_rna_div.py          dxy, pi and net divergence, B and host, all exon positions  -> hostall/01_rna_div_exon.txt
#  7.2 b_rna_div.py (FF=1)   the same at four-fold sites                                -> hostall/02_rna_div_ff.txt
#  7.3 b_rna_boot.py         block bootstrap of the B / host ratios                     -> hostall/03_boot_ff.txt, 04_boot_exon.txt
#  7.4 b_tree_dist.py, b_tanglegram.R, b_cross_test.R   B and host trees of the carriers, Mantel tests -> hostall/05-07
#  7.5 div_matched.sh -> b_div_matched.py   matched comparison: B genes vs their own donor genes -> hostall/10_div_matched.txt
#  7.6 EQUAL=1               each carrier weighted equally                              -> eqweight/
#  7.7 LEAK                  sensitivity to donor sites with possible B reads           -> leak_sensitivity.txt
#  7.8 loo_carrier.sh        leave-one-carrier-out (13 runs + the full set)             -> loo/loo_summary.tsv
source "$(dirname "$0")/../config.sh"; set -euo pipefail
need samtools bcftools python3 Rscript
D=rna_var; mkdir -p hostall eqweight
step() { local n=$1; shift; echo "== $n"; "$@" > hostall/$n.txt 2>&1 || { echo "   FAILED (see hostall/$n.txt)"; return 1; }; }

step 01_rna_div_exon python3 $S/b_rna_div.py $D
step 02_rna_div_ff   env FF=1 python3 $S/b_rna_div.py $D
step 03_boot_ff      python3 $S/b_rna_boot.py $D ff
step 04_boot_exon    python3 $S/b_rna_boot.py $D exon
grep -h "B / host  RNA" hostall/01_rna_div_exon.txt hostall/02_rna_div_ff.txt || true
step 05_tree_dist    python3 $S/b_tree_dist.py $D
step 06_tanglegram   Rscript $S/b_tanglegram.R
step 07_cross_test   Rscript $S/b_cross_test.R
step 10_div_matched  bash $S/div_matched.sh
grep -E "CNsea vs EUinv" hostall/10_div_matched.txt | tail -3 || true

echo "== 7.6 equal weighting"
env FF=1 EQUAL=1 python3 $S/b_rna_div.py $D > eqweight/ff_eq1_7CNsea.txt 2>&1
env EQUAL=1 python3 $S/b_rna_boot.py $D ff > eqweight/boot_ff_eq1.txt 2>&1
env EQUAL=1 python3 $S/b_div_matched.py stats matched > eqweight/matched_eq1.txt 2>&1
grep -h "B / host  RNA" eqweight/ff_eq1_7CNsea.txt || true

echo "== 7.7 donor sites with possible B reads"
for setting in "none 0" "0.1 0" "0.2 0" "0.3 0" "0.2 1"; do
  set -- $setting
  echo "########## LEAK=$1 both directions=$2"
  LEAK=$1 LEAK_BOTHWAYS=$2 python3 $S/b_div_matched.py stats matched 2>&1 | grep -E "sites removed|ratio|CNsea.*EUinv|EUinv.*Med|CNother" || true
done > leak_sensitivity.txt || true

echo "== 7.8 leave-one-carrier-out"
bash $S/loo_carrier.sh
echo "=== step 7 done"
