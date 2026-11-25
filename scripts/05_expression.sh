#!/bin/bash
# Step 5. Co-expression and trans effects (Methods: "Expression, co-expression and trans effects"; Figure 4, Table S5)
#  5.1 wgcna_tetra_excl_gf.R   WGCNA on genotype means (47 genotypes, NAnat excluded) -> wgcna_excl_gf/
#  5.2 b_module_carrier.R      module eigengene ~ lineage + carrier -> wgcna_excl_gf/module_carrier.tsv
#  5.3 bchrom_v4_noPJ.R        carrier vs non-carrier DE (pseudo-bulk, ~ lineage + carrier), trans candidates
#  5.4 bchrom_v4_noPJ.R x 6, trans_lolo_summary.R   leave-one-lineage-out stability of the trans candidates
source "$(dirname "$0")/../config.sh"; set -euo pipefail
mkdir -p logs
echo "== 5.1"; Rscript $S/wgcna_tetra_excl_gf.R > logs/wgcna.log 2>&1; grep -iE "genotypes|soft|modules" logs/wgcna.log | head || true
echo "== 5.2"; Rscript $S/b_module_carrier.R
echo "== 5.3"; Rscript $S/bchrom_v4_noPJ.R > logs/bchrom_v4_noPJ.log 2>&1; grep -E "carriers:|trans candidates" logs/bchrom_v4_noPJ.log | head || true
echo "== 5.4"
for p in CNsea CNfresh CNlake Med EU EUinv; do
  EXCLUDE_POPS=$p OUTDIR=lolo_$p Rscript $S/bchrom_v4_noPJ.R > logs/lolo_$p.log 2>&1 || echo "run without $p failed (logs/lolo_$p.log)"
done
Rscript $S/trans_lolo_summary.R
echo "=== step 5 done"
