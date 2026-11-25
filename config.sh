# Paths and settings for all steps. Edit this file, then run the steps in order (see README.md).
# Every step script reads it with:  source "$(dirname "$0")/../config.sh"

## project folder: one sub-folder per RNA library (see README, "Input layout"); all outputs are written here
export PROJECT=/scratch/project_2009273/Cui/Salt/SaltRNA/Tetra

## reference genome and annotation (Wang et al., 2024; subgenome-resolved assembly that includes PaChr24B)
export GENOME=/scratch//Cui/Salt/annotation_subgenome/genome.subgenome.fa                    # indexed with samtools faidx
export GFF=/scratch/Cui/Salt/annotation_subgenome/07_merge/reed.exon_gene.gff3
export FUNANNO=/scratch/Cui/Salt/annotation_subgenome/07_merge/reed.functional_annotation.tsv
export GOFILE=/scratch/Cui/Salt/annotation_subgenome/07_merge/reed.GO.term2gene.tsv

## nine published genomes: one folder per sample with <sample>/*_1.fastq.gz and <sample>/*_2.fastq.gz
export WGS9_DIR=/scratch/WGS
export WGS9_SAMPLES="EU60GB EU67 EU78 EU620 NAint61 NAint113 NAint186 NAint191 NAint199"

## host exon regions handled separately in variant calling (files in data/; see README)
export LOWDEPTH_BED=${LOWDEPTH_BED:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/data/lowdepth_host.bed}   # called with -d 250
export EXCLUDE_BED=${EXCLUDE_BED:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/data/excluded_host.bed}    # left out

## resources
export THREADS=${THREADS:-8}          # threads per program
export PAR=${PAR:-4}                  # parallel jobs (variant calling per chromosome; leave-one-out runs need ~20 GB each)
export R_LIBS_EXTRA=${R_LIBS_EXTRA:-} # optional extra R library folder

## do not edit below
export REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export S=$REPO/scripts
export BCHR=PaChr24B
need() { for t in "$@"; do command -v $t >/dev/null || { echo "missing program: $t"; exit 1; }; done; }
cd "$PROJECT"
