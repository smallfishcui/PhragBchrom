#!/bin/bash
# Step 1. B transcripts and carrier calls (Methods: "B carrier status", "B transcripts"; Table S1, Table S4)
#  1.1 b_groups.R                genotype / library tables; provisional carriers (% of counts on PaChr24B genes)
#  1.2 assemble                  StringTie on PaChr24B in carriers, merge, gffcompare, quantify in all genotypes
#                                -> stringtie_B_summary.py (B-specific transcripts, ORFs, classes)
#  1.3 rna_carrier_twostep.py    two-step RNA carrier calls -> B_carriers_RNA.tsv
#  1.4 b_groups.R                tables again with the two-step calls (18 carriers)
#  1.5 assemble                  again in the 18 carriers; all later steps use this assembly
#  1.6 salt libraries            B transcripts quantified in the salt library of every genotype
#  1.7 b_q18_check.py            source of the excess B reads in the Q18 control library
#  (carrier_stability.py, which tabulates the carrier calls per library, needs the genome results and runs in step 11)
source "$(dirname "$0")/../config.sh"; set -euo pipefail
need samtools stringtie gffcompare python3 Rscript
O=stringtie_B; mkdir -p $O/bam $O/asm $O/quant

## B sequence and its reference annotation (GFF3 exons -> GTF)
samtools faidx $GENOME $BCHR > $O/B.fa; samtools faidx $O/B.fa
awk -F'\t' -v b=$BCHR '$1==b' $GFF > $O/ref_B.gff3
awk -F'\t' -v OFS='\t' '$3=="exon"{
    t=$9; sub(/.*Parent=/,"",t); sub(/;.*/,"",t)
    g=$9; if (g ~ /gene_id=/) {sub(/.*gene_id=/,"",g); sub(/;.*/,"",g)} else {g=t; sub(/\.[0-9]+$/,"",g)}
    print $1,"ref","exon",$4,$5,".",$7,".","gene_id \""g"\"; transcript_id \""t"\";"
    if(!(t in s) || $4<s[t]) s[t]=$4; if($5>e[t]) e[t]=$5; c[t]=$1; st[t]=$7; gg[t]=g }
  END{ for(t in s) print c[t],"ref","transcript",s[t],e[t],".",st[t],".","gene_id \""gg[t]"\"; transcript_id \""t"\";" }' \
  $O/ref_B.gff3 | sort -k4,4n -k3,3r > $O/ref_B.gtf

assemble() {
  # carriers / non-carriers from column 8 (carrier) of B_groups_genotype.tsv; control library = genotype name
  awk -F'\t' 'NR>1 && $8=="TRUE"  {print $1"\t"$2"\tNA"}' B_groups_genotype.tsv > $O/carriers.txt
  awk -F'\t' 'NR>1 && $8=="FALSE" {print $1"\t"$2"\t0"}'  B_groups_genotype.tsv > $O/noncarriers.txt
  cat $O/carriers.txt $O/noncarriers.txt | cut -f1 > $O/all_samples.txt
  echo "carriers: $(wc -l < $O/carriers.txt)  non-carriers: $(wc -l < $O/noncarriers.txt)"
  rm -f $O/asm/*.gtf $O/quant/*.gtf
  while read -r s; do
    samtools view -b -q 255 -F 0x904 $s/Aligned.sorted.bam $BCHR > $O/bam/$s.bam
    echo -e "$s\t$(grep 'Uniquely mapped reads number' $s/Log.final.out | cut -f2)"
  done < $O/all_samples.txt > $O/libsize.txt
  cut -f1 $O/carriers.txt | while read -r s; do
    stringtie $O/bam/$s.bam -G $O/ref_B.gtf -p $THREADS -m 200 -c 2.5 -f 0.1 -j 2 -l $s -o $O/asm/$s.gtf
  done
  ls $O/asm/*.gtf > $O/asm_list.txt
  stringtie --merge -G $O/ref_B.gtf -m 200 -c 2 -F 1 -T 1 -f 0.05 -o $O/merged_raw.gtf $O/asm_list.txt
  # drop merged transcripts with an intron > 50 kb (spurious long-range joins, e.g. between the two SCC3 copies)
  awk -F'\t' -v M=50000 '$3=="exon"{match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART,RLENGTH);
       if((t in last) && $4-last[t]-1>M) bad[t]=1; last[t]=$5} END{for(t in bad) print t}' $O/merged_raw.gtf > $O/long_intron_tx.txt
  grep -vFf $O/long_intron_tx.txt $O/merged_raw.gtf > $O/merged.gtf || cp $O/merged_raw.gtf $O/merged.gtf
  gffcompare -r $O/ref_B.gtf -o $O/cmp $O/merged.gtf
  while read -r s; do stringtie $O/bam/$s.bam -e -G $O/merged.gtf -p $THREADS -o $O/quant/$s.gtf; done < $O/all_samples.txt
  python3 $S/stringtie_B_summary.py $O
  rm -f $O/bam/*.bam
}

echo "== 1.1"; rm -f B_carriers_RNA.tsv; Rscript $S/b_groups.R
echo "== 1.2"; assemble
echo "== 1.3"; python3 $S/rna_carrier_twostep.py $O
echo "== 1.4"; Rscript $S/b_groups.R
echo "== 1.5"; assemble

echo "== 1.6"   # salt library = genotype name + "S"
mkdir -p $O/bam_salt $O/quant_salt; : > $O/libsize_salt.txt
awk -F'\t' 'NR>1{print $1}' B_groups_genotype.tsv | while read -r g; do
  s=${g}S; [ -s $s/Aligned.sorted.bam ] || { echo "no salt library for $g"; continue; }
  samtools view -b -q 255 -F 0x904 $s/Aligned.sorted.bam $BCHR > $O/bam_salt/$g.bam
  echo -e "$g\t$(grep 'Uniquely mapped reads number' $s/Log.final.out | cut -f2)" >> $O/libsize_salt.txt
  stringtie $O/bam_salt/$g.bam -e -G $O/merged.gtf -p $THREADS -o $O/quant_salt/$g.gtf
done
rm -rf $O/bam_salt

echo "== 1.7"; python3 $S/b_q18_check.py $O Q18 61 97
echo "=== step 1 done"
