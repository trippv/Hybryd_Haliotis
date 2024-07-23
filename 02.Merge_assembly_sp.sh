#!/bin/bash
#SBATCH --job-name=Merge-Assembly
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rgomez@uabc.edu.mx

CPU=$SLURM_NPROCS

REF_GTF_FILE=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Haliotis/Hybryd/RF_Ref.fna

EXPORT=/LUSTRE/apps/bioinformatica/stringtie/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/gffread
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/subread-2.0.6-Linux-x86_64/bin
export PATH=$PATH:$EXPORT

WD=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Haliotis/05.Meta_assembly

ls -d -1 $WD/*_cross/S2_STRINGTIE_REFBASED_MODE/*_transcripts.gtf > stringtie_gtf_list.txt

if [ ! -f "transcripts.gtf" ]; then
 stringtie --rf --merge -p $CPU -o transcripts.gtf stringtie_gtf_list.txt
 
 gffread -w transcripts.fa -g $REF_GTF_FILE transcripts.gtf

else
    echo "Merging step already exists. Continue w/ Feature Count"
fi

# stringtie --rf --merge -p $CPU -o transcripts.gtf stringtie_gtf_list.txt

# gffread -w transcripts.fa -g $REF_GTF_FILE transcripts.gtf

# 2) Count 
$WD/*_cross/S1_HISAT2_SAM_BAM_FILES/*sorted.bam

ls -d -1 $WD/*_cross/S1_HISAT2_SAM_BAM_FILES/*sorted.bam > sorted_bam_list.txt

featureCounts -T $CPU -a transcripts.gtf -o feature_counts.txt sorted_bam_list.txt

exit
