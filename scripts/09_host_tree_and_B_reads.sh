#!/bin/bash
# Step 9. Host tree of all genotypes and distance of the B to host lineages (Methods: "Host tree of all genotypes",
#         "Distance of the B to host lineages"; Figure 1a, Table S1, Table S7). Needs step 6.
#  9.1 host_tree51.py           host distances among all 51 genotypes from RNA SNPs, 100 block bootstraps (needs numpy;
#                               ~48 GB memory) -> hostall/08_host_tree51.txt, figures/data/host51_*.tsv
#  9.2 host_tree_check.R        NJ tree, monophyly and support of each lineage -> hostall/09_host_tree.txt
#  9.3 b_pj_test.py             reads on the B and mismatch rates per library; NAnat-specific loci -> supp_tables/S6_pj_test.txt
#  9.4 b_spec_lineage_check.py  B-specific transcripts expressed by non-carriers of one continent
source "$(dirname "$0")/../config.sh"; set -euo pipefail
need samtools python3 Rscript
mkdir -p hostall figures/data supp_tables logs
echo "== 9.1"; python3 $S/host_tree51.py rna_var 100 > hostall/08_host_tree51.txt 2>&1; tail -3 hostall/08_host_tree51.txt
echo "== 9.2"; Rscript $S/host_tree_check.R > hostall/09_host_tree.txt 2>&1
echo "== 9.3"; python3 $S/b_pj_test.py stringtie_B > supp_tables/S6_pj_test.txt 2>&1
echo "== 9.4"; python3 $S/b_spec_lineage_check.py stringtie_B > logs/spec_lineage.txt 2>&1
echo "=== step 9 done"
