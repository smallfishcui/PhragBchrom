# Small input files

|File|What it is|Used in|
|-|-|-|
|`lowdepth\_host.bed`|Host exon intervals (BED: chrom, start, end) where variants are called with at most 250 reads per library and position (`bcftools mpileup -d 250`) instead of 5,000. These are very deep regions where calling with `-d 5000` does not finish.|step 6.3|
|`excluded\_host.bed`|Two PaChr5A regions (66 exon intervals, 36 kb, 0.07% of host exon sequence) left out of both the host variants and the callable positions.|steps 6, 7|

Both paths are set in `config.sh` (`LOWDEPTH\_BED`, `EXCLUDE\_BED`). Leave a file out to call all regions with the
default settings; results then differ slightly from the paper.

