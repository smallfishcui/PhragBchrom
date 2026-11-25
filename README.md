# The B chromosome of the common reed (*Phragmites australis*): analysis code

Code for the analyses, figures and supplementary tables of
"*Diversification and spread of B chromosomes across host lineages of the common reed*". All analyses use RNA-seq of 51 genotypes (one control and one
salt-treated library each), the published reference genome and annotation, and nine published genomes.

## Layout

```
config.sh        paths and settings (edit this first)
steps/           one script per analysis step, run in order 01 -> 11
scripts/         the Python and R scripts the steps call (each script names its step in its first lines)
data/            small input files that are not produced by the code (see data/README.md)
```

## How to run

1. Install the software below.
2. Edit `config.sh`: project folder, reference files, folder with the nine published genomes, threads.
3. Put the inputs in the project folder (layout below).
4. Run the steps in order. Each step writes into the project folder and reads what earlier steps wrote:

```bash
bash steps/01\\\_B\\\_transcripts\\\_and\\\_carriers.sh
bash steps/02\\\_published\\\_genomes.sh
...
bash steps/11\\\_figures\\\_tables.sh
```

Steps 2, 3, 4 and 5 depend only on step 1 and can run in any order. Step 6 (variant calling on all host
chromosomes) is the slowest: hours to days depending on `PAR` and `THREADS`.
The 1,000 salt permutations in step 10.2 can also be split into separate jobs (`PERM\\\_SEED` 1-100, `NPERM=10` each).

## Steps

|Step|Script in `steps/`|Sub-step: script in `scripts/`|Manuscript|
|-|-|-|-|
|1|`01\\\_B\\\_transcripts\\\_and\\\_carriers.sh`|1.1 `b\\\_groups.R` (provisional carriers) · 1.2 StringTie assembly on PaChr24B → `stringtie\\\_B\\\_summary.py` · 1.3 `rna\\\_carrier\\\_twostep.py` (two-step carrier calls) · 1.4 `b\\\_groups.R` (final carriers) · 1.5 assembly in the 18 carriers → `stringtie\\\_B\\\_summary.py` · 1.6 quantification in the salt libraries · 1.7 `b\\\_q18\\\_check.py`|B carrier status; B transcripts; Table S1, S4|
|2|`02\\\_published\\\_genomes.sh`|2.1 bwa index · 2.2 mapping, duplicate marking, 1-kb depth · 2.3 mappability of PaChr24B · 2.4 `wgs9\\\_B.py` (B copy number, structure, SKA3 insertion) · 2.5 `b\\\_seg\\\_rna.py`|B copy number and structure; Fig. 1c, 2, 6c|
|3|`03\\\_gene\\\_origin.sh`|3.1 `b\\\_origin.py` (donor genes, Ks, Ka, families) · 3.2 `b\\\_host\\\_control.py` · 3.3 `b\\\_kaks\\\_sens.py` · 3.4 `b\\\_origin\\\_categories.py` · 3.5 `b\\\_keygene\\\_phylo.py`, `codon\\\_thread.py`, `keygene\\\_phylo\\\_summary.R` (gene trees) · 3.6 `b\\\_keygene\\\_structure.py`|Origin and age of B genes; Tables 2–3, S2; Fig. 3a, c, d|
|4|`04\\\_GO.sh`|4.1 `b\\\_go\\\_main.R`|Functional enrichment; Table S3; Fig. 3b|
|5|`05\\\_expression.sh`|5.1 `wgcna\\\_tetra\\\_excl\\\_gf.R` · 5.2 `b\\\_module\\\_carrier.R` · 5.3 `bchrom\\\_v4\\\_noPJ.R` (carrier DE, trans candidates) · 5.4 leave-one-lineage-out → `trans\\\_lolo\\\_summary.R`|Expression, co-expression and trans effects; Fig. 4; Table S5|
|6|`06\\\_variants.sh`|6.1–6.3 bcftools calls on B and host exons · 6.4 focal carriers · 6.5 `b\\\_fourfold.py` · 6.6 callable positions · 6.7 SKA3 insertion reads|B types and divergence (Methods)|
|7|`07\\\_divergence.sh`|7.1–7.2 `b\\\_rna\\\_div.py` · 7.3 `b\\\_rna\\\_boot.py` · 7.4 `b\\\_tree\\\_dist.py`, `b\\\_tanglegram.R`, `b\\\_cross\\\_test.R` · 7.5 `div\\\_matched.sh` → `b\\\_div\\\_matched.py` · 7.6 equal weighting · 7.7 donor-site filter sensitivity · 7.8 `loo\\\_carrier.sh`|B types and divergence; Table 4; Fig. 7; Table S6|
|8|`08\\\_splicing.sh`|8.1 `sj\\\_B.R` · 8.2 `stringtie\\\_B\\\_composition.py` · 8.3 `b\\\_isoform\\\_switch\\\_isa.R` · 8.4 `sj\\\_variants\\\_rna.py` · 8.5 `ska3\\\_table.py` · 8.6 `igv\\\_regions.py`, `sashimi\\\_data.py`, `sashimi\\\_B.R` · 8.7 `b\\\_salt\\\_validate.py`|Splice junctions; splice-site variants; SKA3; Fig. 5, 6; Table S4|
|9|`09\\\_host\\\_tree\\\_and\\\_B\\\_reads.sh`|9.1 `host\\\_tree51.py` · 9.2 `host\\\_tree\\\_check.R` · 9.3 `b\\\_pj\\\_test.py` · 9.4 `b\\\_spec\\\_lineage\\\_check.py`|Host tree; distance of the B to host lineages; Fig. 1a; Tables S1, S7|
|10|`10\\\_salt.sh`|10.1 `b\\\_salt\\\_B.R` · 10.2 `b\\\_salt\\\_B.R` (PERM\_ONLY, 1,000 permutations) · 10.3 `salt\\\_perm\\\_summary.R` · 10.4 `b\\\_salt\\\_Btype.R`|Salt response; Fig. S1; Table S5|
|11|`11\\\_figures\\\_tables.sh`|11.1 `fig\\\_data.py` · 11.2 `carrier\\\_stability.py` · 11.3 `b\\\_mix\\\_check\\\_rna.py` · 11.4 `figures.R` · 11.5 `fig8\\\_model.R` · 11.6 `make\\\_supp\\\_tables.py`|Figures 1–8, S1; Tables S1–S7|

The comment block at the top of each step script lists its sub-steps, inputs and outputs.

## Input layout (project folder)

RNA-seq reads were trimmed with Trimmomatic, aligned with STAR 2.7.11b, sorted with samtools and counted with
htseq-count 2.1.2 (see Methods). These steps are not included here. The code expects one folder per library:

```
<library>/Aligned.sorted.bam (+ .bai)   STAR alignment, coordinate-sorted (unique reads have MAPQ 255)
<library>/Log.final.out                  STAR summary (uniquely mapped reads are read from it)
<library>/SJ.out.tab                     STAR splice junctions
<library>/<library>\\\_counts.txt           htseq-count table
sample\\\_info.txt                          one row per library: columns "sample" (library name) and "condition" (lineage)
id\\\_map.tsv                               published genome ID -> RNA genotype ("wgs\\\_id", "rna\\\_id")
```

The control library of a genotype is named after the genotype (e.g. `41`), the salt-treated library adds `S` (`41S`).
Lineage codes in `sample\\\_info.txt`: CNfresh, CNlake, CNsea, YH, Med, EU, EUinv, and PJ for the North American
native lineage (NAnat in the paper). NAnat is left out of the carrier analyses.

Reference files (set in `config.sh`): genome `genome.subgenome.fa` (with the B chromosome PaChr24B), annotation
`reed.exon\\\_gene.gff3`, `reed.functional\\\_annotation.tsv` and `reed.GO.term2gene.tsv` (Wang et al., 2024).

## Software

bash, awk, samtools, bcftools, StringTie, gffcompare, bwa, mosdepth, BLAST+ (blastp), MAFFT, IQ-TREE (2 or 3);
Python 3 with numpy, pandas and openpyxl; R with ape, phangorn, DESeq2, apeglm, ashr, WGCNA, clusterProfiler, GO.db
and IsoformSwitchAnalyzeR. Versions are given in the Methods.

