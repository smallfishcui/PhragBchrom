#!/bin/bash
# Step 6. SNPs from the RNA reads and callable positions (Methods: "B types and divergence")
#  6.1 BAM list (control and salt library of every genotype; read counts are summed per genotype later)
#  6.2 regions: exons of B-specific transcripts (PaChr24B); exons of host genes on all 24 A/D chromosomes
#  6.3 bcftools mpileup/call per chromosome (-q 255 -Q 20 -I -d 5000). Host regions in LOWDEPTH_BED are called with
#      -d 250 (very deep regions that otherwise do not finish); regions in EXCLUDE_BED are left out of both the
#      variants and the callable positions (36 kb, 0.07% of host exon sequence)
#  6.4 focal carriers (7 CNsea, 6 EUinv) whose reads define callable positions
#  6.5 b_fourfold.py: four-fold degenerate sites of B and host genes
#  6.6 callable positions: >= 10 reads summed over the CNsea carriers and over the EUinv carriers
#      (all exon positions and four-fold sites; B and host)
#  6.7 reads at the SKA3 GG insertion
source "$(dirname "$0")/../config.sh"; set -euo pipefail
need samtools bcftools python3
O=rna_var; mkdir -p $O/chr
[ -s $GENOME.fai ] || samtools faidx $GENOME

echo "== 6.1"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="sample") k=i; next} {print $k"/Aligned.sorted.bam"}' B_groups_sample.tsv \
  | while read -r b; do [ -s "$b" ] && echo "$b"; done > $O/bams.txt
echo "BAMs: $(wc -l < $O/bams.txt)"

echo "== 6.2"
merge() { sort -k1,1 -k2,2n | awk -F'\t' -v OFS='\t' '{ if($1==c && $2<=e){ if($3>e) e=$3 } else { if(c!="") print c,s,e; c=$1; s=$2; e=$3 } } END{ if(c!="") print c,s,e }'; }
minus() {   # intervals of $2 that do not overlap any interval of $1
  if [ -s "$1" ]; then awk -F'\t' 'NR==FNR{c[NR]=$1; s[NR]=$2; e[NR]=$3; n=NR; next} {k=1; for(i=1;i<=n;i++) if($1==c[i] && $3>s[i] && $2<e[i]) {k=0; break} if(k) print}' "$1" "$2"
  else cat "$2"; fi; }
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="B_specific") k=i; next} $k=="True"{print $1}' stringtie_B/B_transcripts.tsv > $O/Bspec_tx.txt
awk -F'\t' -v L=$O/Bspec_tx.txt 'BEGIN{while((getline t < L)>0) k[t]=1}
  $3=="exon"{match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16); if(t in k) print $1"\t"$4-1"\t"$5}' \
  stringtie_B/cmp.annotated.gtf | merge > $O/B.bed
grep -P '^PaChr[0-9]+[AD]\t' $GENOME.fai | cut -f1 > $O/host_chroms.txt
awk -F'\t' -v L=$O/host_chroms.txt 'BEGIN{while((getline c < L)>0) k[c]=1} ($1 in k) && $3=="exon"{print $1"\t"$4-1"\t"$5}' $GFF | merge > $O/host.bed
minus "$EXCLUDE_BED" $O/host.bed > $O/host_called.bed
[ -s "$EXCLUDE_BED" ] && cp $EXCLUDE_BED $O/excluded_host.bed      # read by later steps
echo "B-specific exons: $(awk '{s+=$3-$2} END{print s}' $O/B.bed) bp; host exons: $(awk '{s+=$3-$2} END{print s}' $O/host.bed) bp" \
     "($(awk '{s+=$3-$2} END{print s}' $O/host_called.bed) bp called) on $(wc -l < $O/host_chroms.txt) chromosomes"

echo "== 6.3"
call() {   # chromosome, region BED
  local c=$1 reg=$2 d=$O/chr/$1
  rm -f $d.all.bed.*.bed $d.*.bcf
  awk -v c=$c '$1==c' $reg > $d.all.bed
  if [ -s "$LOWDEPTH_BED" ]; then
    awk -F'\t' 'NR==FNR{c[NR]=$1; s[NR]=$2; e[NR]=$3; n=NR; next} {k=0; for(i=1;i<=n;i++) if($1==c[i] && $3>s[i] && $2<e[i]) {k=1; break} print > (FILENAME"." (k ? "low" : "std") ".bed")}' \
        $LOWDEPTH_BED $d.all.bed
  else cp $d.all.bed $d.all.bed.std.bed; fi
  : > $d.ad.tsv
  for part in std low; do
    f=$d.all.bed.$part.bed; [ -s $f ] || continue
    bcftools mpileup -f $GENOME -R $f -b $O/bams.txt --ignore-RG -q 255 -Q 20 -I -d $([ $part = low ] && echo 250 || echo 5000) \
        -a FORMAT/AD,FORMAT/DP -Ou 2>/dev/null | bcftools call -mv -Ob -o $d.$part.bcf 2>/dev/null
    bcftools query -f '%POS\t%REF\t%ALT[\t%AD]\n' $d.$part.bcf >> $d.ad.tsv
  done
  sort -k1,1n -u $d.ad.tsv > $O/$c.ad.tsv
  echo "$c: $(wc -l < $O/$c.ad.tsv) variant sites"
}
export -f call; export O GENOME LOWDEPTH_BED
( echo "$BCHR $O/B.bed"; while read c; do echo "$c $O/host_called.bed"; done < $O/host_chroms.txt ) | xargs -P $PAR -L 1 bash -c 'call $0 $1'
bcftools query -l $O/chr/$BCHR.std.bcf > $O/samples_order.txt

echo "== 6.4"
awk -F'\t' 'NR>1 && $8=="TRUE" && ($2=="CNsea"||$2=="EUinv"){print $1"\t"$2}' B_groups_genotype.tsv > $O/focal.txt
python3 - "$O" <<'PY'
import sys; D = sys.argv[1]
foc = [l.split() for l in open(f"{D}/focal.txt")]
with open(f"{D}/focal_bams.txt", "w") as o, open(f"{D}/focal_groups.txt", "w") as g:
    for b in (l.strip() for l in open(f"{D}/bams.txt")):
        s = b.split("/")[0]; geno = s[:-1] if s.endswith("S") else s
        for f, p in foc:
            if f == geno: o.write(b + "\n"); g.write(p + "\n")
PY
echo "focal carriers: $(cut -f2 $O/focal.txt | sort | uniq -c | tr '\n' ' '); libraries: $(wc -l < $O/focal_bams.txt)"

echo "== 6.5"
python3 $S/b_fourfold.py $GENOME $GFF $O
minus "$EXCLUDE_BED" $O/ff_host.bed > $O/ff_host.tmp && mv $O/ff_host.tmp $O/ff_host.bed

echo "== 6.6"
GRP=$(paste -sd, $O/focal_groups.txt)
for job in "exon B $O/B.bed" "exon host $O/host_called.bed" "ff B $O/ff_B.bed" "ff host $O/ff_host.bed"; do
  set -- $job
  samtools depth -a -Q 255 -b $3 -f $O/focal_bams.txt 2>/dev/null \
  | awk -v G=$GRP 'BEGIN{n=split(G,g,",")} {cn=0; eu=0; for(i=1;i<=n;i++){ if(g[i]=="CNsea") cn+=$(i+2); else eu+=$(i+2) } if(cn>=10 && eu>=10) print $1"\t"$2}' \
  > $O/callable_pos_$1_$2.txt
  if [ $1 = ff ]; then f=callable_ff_$2.txt; else f=callable_$2.txt; fi
  wc -l < $O/callable_pos_$1_$2.txt > $O/$f
  echo "callable $1 $2: $(cat $O/$f)"
done

echo "== 6.7"
samtools mpileup -r $BCHR:20171600-20171700 -q 255 -Q 0 -B -b $O/bams.txt 2>/dev/null > $O/ska3.mpileup
echo "=== step 6 done"
