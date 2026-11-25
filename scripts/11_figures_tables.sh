#!/bin/bash
# Step 11. Figures 1-8, S1 and Supplementary Tables S1-S7. Needs steps 1-10.
#  11.1 fig_data.py            figure data tables -> figures/data/
#  11.2 carrier_stability.py   carrier calls per library, agreement with the published genomes -> carrier_stability.tsv
#  11.3 b_mix_check_rna.py     mixed B types within a carrier (diagnostic sites) -> figures/data/mix_rna.txt
#  11.4 figures.R              Figures 1-7 and S1 -> figures/
#  11.5 fig8_model.R           Figure 8 (schematic)
#  11.6 make_supp_tables.py    Tables S1-S7 (one Excel file, one sheet per source table) -> supp_tables/
source "$(dirname "$0")/../config.sh"; set -euo pipefail
mkdir -p figures/data supp_tables
nb=$(cut -f1 figures/data/host51_boot.tsv 2>/dev/null | sort -u | grep -vc rep || true)
[ "${nb:-0}" -ge 100 ] || { echo "figures/data/host51_boot.tsv has ${nb:-0} bootstrap replicates: run step 9 first"; exit 1; }
echo "== 11.1"; python3 $S/fig_data.py
echo "== 11.2"; python3 $S/carrier_stability.py stringtie_B
echo "== 11.3"; python3 $S/b_mix_check_rna.py rna_var > figures/data/mix_rna.txt
echo "== 11.4"; Rscript $S/figures.R
echo "== 11.5"; Rscript $S/fig8_model.R figures
echo "== 11.6"; python3 $S/make_supp_tables.py
echo "=== step 11 done"
