#!/bin/bash
#SBATCH --job-name=hisat2-build
#SBATCH -N 1
#SBATCH --mem=80GB
#SBATCH --ntasks-per-node=20

CPU=$SLURM_NPROCS

genome=$1 # 
prefix=`basename ${genome%.f}`

hisat2-build -p $CPU $genome $prefix

exit