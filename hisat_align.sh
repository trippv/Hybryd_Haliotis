#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH -N 1
#SBATCH --mem=40GB
#SBATCH --ntasks-per-node=20

idx_file=$1
FirstPair=$2 # CFF1_R1_p.fq.gz

CPU=$SLURM_NPROCS

mkdir -p HISAT2_SAM_BAM_FILES

base="${FirstPair%*_*R1_p.fq.gz}"

hisat2  --phred33 --dta -p $CPU \
        -x $idx_file -1 ${base}_R1_p.fq.gz -2 ${base}_R2_p.fq.gz \
        --rg-id=${base} --rg SM:${base} -S HISAT2_SAM_BAM_FILES/${base}.sam \
        --summary-file ${base}.summary.txt --met-file ${base}.met.txt

samtools sort -@ $CPU -o HISAT2_SAM_BAM_FILES/${base}.sorted.bam HISAT2_SAM_BAM_FILES/${base}.sam
