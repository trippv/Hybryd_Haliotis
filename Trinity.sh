#!/bin/bash
#SBATCH --job-name=Trinity
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

module load trinityrnaseq-v2.15.1

Trinity --seqType fq --max_memory 100G --left CFF1_R1_p.fq.gz,CFF2_R1_p.fq.gz,CFF3_R1_p.fq.gz,CRF2_R1_p.fq.gz,CRF3_R1_p.fq.gz,CRR1_R1_p.fq.gz,CRR2_R1_p.fq.gz,CRR3_R1_p.fq.gz --right CFF1_R2_p.fq.gz,CFF2_R2_p.fq.gz,CFF3_R2_p.fq.gz,CRF2_R2_p.fq.gz,CRF3_R2_p.fq.gz,CRR1_R2_p.fq.gz,CRR2_R2_p.fq.gz,CRR3_R2_p.fq.gz --CPU 20 -min_kmer_cov 2


# min_kmer_cov :min count for K-mers to be assembled by Inchworm (default: 1)

exit