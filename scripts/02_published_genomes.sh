#!/bin/bash
# Step 2. Nine published genomes: B copy number and structure (Methods: "B copy number and structure";
#         Figure 1c, Figure 2, Figure 6c, Table S1)
#  2.1 bwa index of the reference (skipped if an index exists next to the genome)
#  2.2 per genome: bwa mem -> fixmate -> sort -> markdup; 1-kb depth at MAPQ >= 20 (mosdepth); small BAM of PaChr24B
#  2.3 mappability of PaChr24B: 150-bp reads every 50 bp mapped back to the genome
#  2.4 SKA3 insertion reads; wgs9_B.py: relative coverage, B copy number, per-copy coverage, absent segments
#  2.5 b_seg_rna.py: RNA check of the segments missing from the European-type B
source "$(dirname "$0")/../config.sh"; set -euo pipefail
need bwa samtools mosdepth python3
W=wgs9; mkdir -p $W/ref $W/bam $W/depth $W/tmp

echo "== 2.1"
ln -sf $GENOME $W/ref/genome.fa
if [ -s $GENOME.bwt ]; then for x in amb ann bwt pac sa; do ln -sf $GENOME.$x $W/ref/genome.fa.$x; done
else bwa index $W/ref/genome.fa; fi
[ -s $W/ref/genome.fa.fai ] || samtools faidx $W/ref/genome.fa
REF=$W/ref/genome.fa

echo "== 2.2"
for s in $WGS9_SAMPLES; do
  [ -s $W/bam/$s.B.bam ] && { echo "$s done"; continue; }
  r1=$(ls $WGS9_DIR/$s/*_1.fastq.gz | head -1); r2=$(ls $WGS9_DIR/$s/*_2.fastq.gz | head -1)
  bwa mem -t $THREADS -R "@RG\tID:$s\tSM:$s" $REF $r1 $r2 2> $W/tmp/$s.bwa.log \
    | samtools fixmate -m -u - - | samtools sort -u -@ 4 -m 2G -T $W/tmp/$s - \
    | samtools markdup -@ 4 - $W/bam/$s.full.bam
  samtools index $W/bam/$s.full.bam
  samtools flagstat $W/bam/$s.full.bam > $W/bam/$s.flagstat
  mosdepth -t 4 -Q 20 -n -F 1796 --by 1000 $W/depth/$s $W/bam/$s.full.bam
  samtools view -b $W/bam/$s.full.bam $BCHR > $W/bam/$s.B.bam; samtools index $W/bam/$s.B.bam
  [ "${KEEP:-0}" = 1 ] || rm -f $W/bam/$s.full.bam $W/bam/$s.full.bam.bai
done

echo "== 2.3"
samtools faidx $REF $BCHR > $W/tmp/B.fa
python3 - <<'PY'
seq = "".join(l.strip() for l in open("wgs9/tmp/B.fa") if not l.startswith(">")).upper()
with open("wgs9/tmp/B_tiles.fq", "w") as o:
    for i in range(0, len(seq) - 150 + 1, 50):
        r = seq[i:i + 150]
        if r.count("N") > 15: continue
        o.write(f"@t{i}\n{r}\n+\n{'I' * 150}\n")
PY
bwa mem -t $THREADS $REF $W/tmp/B_tiles.fq 2>/dev/null | samtools sort -@ 4 -o $W/tmp/B_tiles.bam
samtools index $W/tmp/B_tiles.bam
mosdepth -t 4 -Q 20 -n --chrom $BCHR --by 1000 $W/mappability $W/tmp/B_tiles.bam

echo "== 2.4"
ls $W/bam/*.B.bam > $W/ska3_bams.txt
sed 's#.*/##; s#\.B\.bam$##' $W/ska3_bams.txt > $W/ska3_order.txt
samtools mpileup -r $BCHR:20171600-20171700 -q 20 -Q 0 -B -b $W/ska3_bams.txt 2>/dev/null > $W/ska3.mpileup
python3 $S/wgs9_B.py $W

echo "== 2.5"; python3 $S/b_seg_rna.py stringtie_B
echo "=== step 2 done"
