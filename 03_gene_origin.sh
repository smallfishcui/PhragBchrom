#!/bin/bash
# Step 3. Origin and age of B genes (Methods: "Origin and age of B genes"; Tables 2-3, Table S2, Figure 3a, 3c, 3d)
#  3.1 b_origin.py prep / blastp / analyze   donor (closest host gene) of each B gene, Ks and Ka, donor segments,
#                                            B gene families, A-D homeolog Ks (3,000 A genes sampled)
#  3.2 b_host_control.py                     Ka/Ks controls (host pairs of matched Ks; donor vs its own homeolog)
#  3.3 b_kaks_sens.py                        Ka/Ks with complete, full-length genes only
#  3.4 b_origin_categories.py                Ks and Ka/Ks by functional group (Table 3)
#  3.5 b_keygene_phylo.py, codon_thread.py,  gene trees (MAFFT, codon alignment, IQ-TREE, 1,000 ultrafast bootstraps) of
#      keygene_phylo_summary.R               chromosome-segregation families and HAP2/GEX1-3; placement of the B copies
#  3.6 b_keygene_structure.py                gene structure of CENH3, SCC3, HEI10 and SKA3 copies vs their donors
source "$(dirname "$0")/../config.sh"; set -euo pipefail
need blastp makeblastdb mafft python3 Rscript
IQ=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree || true); [ -n "$IQ" ] || { echo "missing program: iqtree"; exit 1; }
O=B_origin; mkdir -p $O Bgenes
FMT="6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq"

echo "== 3.1"
python3 $S/b_origin.py prep $GENOME $GFF $O $BCHR
awk -F'\t' 'NR>1 && $2 ~ /A$/ {print $1}' $O/genes.tsv | shuf --random-source=<(yes) -n 3000 > $O/A_sample.txt
awk -F'\t' 'NR>1 && $2 ~ /D$/ {print $1}' $O/genes.tsv > $O/D_genes.txt
pick() { awk -v L="$1" 'BEGIN{while((getline g < L)>0) k[g]=1} /^>/{p=(substr($1,2) in k)} p' $O/all.pep.fa; }
pick $O/A_sample.txt > $O/A_sample.pep.fa; pick $O/D_genes.txt > $O/D.pep.fa
makeblastdb -in $O/all.pep.fa -dbtype prot -out $O/allDB >/dev/null
makeblastdb -in $O/D.pep.fa   -dbtype prot -out $O/DDB   >/dev/null
blastp -query $O/B.pep.fa        -db $O/allDB -evalue 1e-10 -max_target_seqs 20 -outfmt "$FMT" -num_threads $THREADS > $O/B_vs_all.tsv
blastp -query $O/A_sample.pep.fa -db $O/DDB   -evalue 1e-10 -max_target_seqs 5  -outfmt "$FMT" -num_threads $THREADS > $O/A_vs_D.tsv
python3 $S/b_origin.py analyze $O $BCHR

echo "== 3.2"
cut -f3 $O/B_origin_pairs.tsv | tail -n +2 | grep . | sort -u > $O/host_genes.txt
pick $O/host_genes.txt > $O/host.pep.fa
blastp -query $O/host.pep.fa -db $O/allDB -evalue 1e-10 -max_target_seqs 20 -outfmt "$FMT" -num_threads $THREADS > $O/host_vs_all.tsv
python3 $S/b_host_control.py $O

echo "== 3.3"; python3 $S/b_kaks_sens.py $O
echo "== 3.4"; python3 $S/b_origin_categories.py $O

echo "== 3.5"
K=keygene_phylo_go; rm -rf $K; mkdir -p $K
python3 $S/b_keygene_phylo.py $O $K
for f in $K/*.pep.fa; do
  s=${f%.pep.fa}; [ $(grep -c ">" $f) -ge 4 ] || { echo "skip $(basename $s): < 4 sequences"; continue; }
  mafft --auto --quiet $f > $s.aln.pep.fa
  python3 $S/codon_thread.py $s.aln.pep.fa $s.cds.fa > $s.aln.cds.fa
  $IQ -s $s.aln.cds.fa -m MFP -B 1000 -T 4 --prefix $s -redo -quiet || echo "iqtree failed: $s"
done
KP_OUT=$K Rscript $S/keygene_phylo_summary.R

echo "== 3.6"; python3 $S/b_keygene_structure.py $GENOME $GFF $FUNANNO > Bgenes/keygene_structure.txt; cat Bgenes/keygene_structure.txt
echo "=== step 3 done"
