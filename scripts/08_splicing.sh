#!/bin/bash
# Step 8. Splice junctions and B-type differences in transcripts (Methods: "Splice junctions", "Sequence variants at
#         splice sites", "SKA3 junctions"; Figure 5, Figure 6, Table S4). Needs steps 1, 2 and 6.
#  8.1 sj_B.R                      B-specific junctions (STAR SJ.out.tab) with >= 2, 3, 5, 10 reads -> sj_B_min<k>/
#                                  (the paper uses sj_B_min3/)
#  8.2 stringtie_B_composition.py  share of B expression per transcript, CNsea vs EUinv (permutation tests)
#  8.3 b_isoform_switch_isa.R      isoform switches CNsea vs EUinv (IsoformSwitchAnalyzeR), control and salt libraries
#  8.4 sj_variants_rna.py          differing SNPs near splice sites, switching vs non-switching loci -> sjvar_isa.txt
#  8.5 ska3_table.py               SKA3 junction reads per library -> ska3_junctions.tsv
#  8.6 igv_regions.py, sashimi_data.py, sashimi_B.R   group BAMs, IGV session and sashimi plot of SKA3 -> igv_B/
#  8.7 b_salt_validate.py          8.2, 8.3 and SKA3 repeated in the salt-treated libraries
source "$(dirname "$0")/../config.sh"; set -euo pipefail
need samtools python3 Rscript
mkdir -p logs

echo "== 8.1"
for k in 2 3 5 10; do MINREADS=$k SJOUT=sj_B_min$k Rscript $S/sj_B.R > logs/sj_B_min$k.log 2>&1 || echo "failed: logs/sj_B_min$k.log"; done
echo "== 8.2"; python3 $S/stringtie_B_composition.py stringtie_B
echo "== 8.3"; READLEN=150 Rscript $S/b_isoform_switch_isa.R
echo "== 8.4"; SWITCH=B_gene_isoform_switch_ISA.tsv MERGE=10 python3 $S/sj_variants_rna.py stringtie_B rna_var > sjvar_isa.txt
echo "== 8.5"; python3 $S/ska3_table.py

echo "== 8.6"
O=igv_B; mkdir -p $O/bam
python3 $S/igv_regions.py stringtie_B $GFF > $O/regions.tsv
REG=$(awk 'NR>1{printf "%s:%s-%s ", $2, $3, $4}' $O/regions.tsv)
grp() { case $1 in
  CNsea|EUinv)    awk -F'\t' -v p=$1 '$2==p{print $1}' stringtie_B/carriers.txt;;
  other_carriers) awk -F'\t' '$2!="CNsea" && $2!="EUinv"{print $1}' stringtie_B/carriers.txt;;
  noncarriers)    cut -f1 stringtie_B/noncarriers.txt;; esac; }
for g in CNsea EUinv other_carriers noncarriers; do
  parts=()
  for s in $(grp $g); do samtools view -b -M -q 255 $s/Aligned.sorted.bam $REG > $O/bam/tmp_$s.bam; parts+=($O/bam/tmp_$s.bam); done
  samtools merge -f -o - "${parts[@]}" | samtools sort -o $O/bam/RNA_$g.bam; samtools index $O/bam/RNA_$g.bam
  rm -f "${parts[@]}"
done
for t in NAint EU; do          # published genomes, SKA3 region
  ls wgs9/bam/${t}*.B.bam > $O/wgs_$t.txt
  samtools merge -f -b $O/wgs_$t.txt -R $BCHR:20160000-20185000 -o - | samtools sort -o $O/bam/WGS_$t.bam; samtools index $O/bam/WGS_$t.bam
done
cp stringtie_B/B.fa stringtie_B/B.fa.fai stringtie_B/ref_B.gtf $O/; cp stringtie_B/merged.gtf $O/stringtie_merged.gtf
{ echo '<?xml version="1.0" encoding="UTF-8"?>'
  echo "<Session genome=\"B.fa\" locus=\"$(awk 'NR==2{print $2":"$3"-"$4}' $O/regions.tsv)\" version=\"8\">"
  echo '  <Resources>'
  for f in RNA_CNsea RNA_EUinv RNA_other_carriers RNA_noncarriers WGS_NAint WGS_EU; do echo "    <Resource path=\"bam/$f.bam\"/>"; done
  echo '    <Resource path="ref_B.gtf"/>'; echo '    <Resource path="stringtie_merged.gtf"/>'
  echo '  </Resources>'; echo '  <Regions>'
  awk 'NR>1{printf "    <RegionOfInterest chromosome=\"%s\" start=\"%d\" end=\"%d\" description=\"%s\"/>\n", $2, $3-1, $4, $1}' $O/regions.tsv
  echo '  </Regions>'; echo '</Session>'; } > $O/B_loci.xml
python3 $S/sashimi_data.py $O
Rscript $S/sashimi_B.R $O

echo "== 8.7"; python3 $S/b_salt_validate.py stringtie_B
echo "=== step 8 done"
