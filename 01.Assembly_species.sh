#!/bin/bash
#SBATCH --job-name=Guide-Assembly
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rgomez@uabc.edu.mx

# To run
# Ex. /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Haliotis/05.Meta_assembly
# sbatch 01.Assembly_species.sh Reference.fasta Annotation.gtf

EXPORT=/LUSTRE/apps/bioinformatica/hisat2/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/stringtie/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/gffread
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/TransDecoder-v5.7.0/
export PATH=$PATH:$EXPORT

# Vars

REFERENCE=$1 # FASTA FILE 

REF_GTF_FILE=$2

REF_PREFIX=`basename ${REFERENCE%.f*}`

CPU=$SLURM_NPROCS


# 1) Create index if does not exist

mkdir -p INDEX

if [ ! -f "INDEX/${REF_PREFIX}.1.ht2" ]; then
hisat2-build -p $CPU $REFERENCE INDEX/$REF_PREFIX

else
    echo "index for '$REFERENCE' already exists."
fi

# 2) Align reads back to the reference

mkdir -p S1_HISAT2_SAM_BAM_FILES

for i in $(ls *_R1_p.fq.gz)
do
bs="${i%*_R1_p.fq.gz}"

# Test if the alignment was previously done!

if [ ! -f "S1_HISAT2_SAM_BAM_FILES/${bs}.sorted.bam" ]; then

left_file=${bs}_R1_p.fq.gz
right_file=${bs}_R2_p.fq.gz

hisat2  --phred33 --dta -p $CPU \
        -x INDEX/$REF_PREFIX -1 $left_file -2 $right_file \
        --rg-id=${bs} --rg SM:${bs} -S S1_HISAT2_SAM_BAM_FILES/${bs}.sam \
        --summary-file ${bs}.summary.txt --met-file ${bs}.met.txt

samtools sort -@ $CPU -o S1_HISAT2_SAM_BAM_FILES/${bs}.sorted.bam S1_HISAT2_SAM_BAM_FILES/${bs}.sam


# 3) Scaffolding reads w/ Stringtie
mkdir -p S2_STRINGTIE_REFBASED_MODE;   

stringtie --rf -p $CPU -G $REF_GTF_FILE -l $bs -o S2_STRINGTIE_REFBASED_MODE/${bs}_transcripts.gtf S1_HISAT2_SAM_BAM_FILES/${bs}.sorted.bam


else
    echo "Alignment for 'S1_HISAT2_SAM_BAM_FILES/${bs}.sorted.bam' already exists."
fi


done

# 4) Get fasta (This is addit. step to comparative transcriptomic step)
ls -1 S2_STRINGTIE_REFBASED_MODE/*_transcripts.gtf > stringtie_gtf_list.tmp

stringtie --rf --merge -p $CPU -o S2_STRINGTIE_REFBASED_MODE/transcripts.gtf stringtie_gtf_list.tmp

gffread -w transcripts.fa -g $REFERENCE S2_STRINGTIE_REFBASED_MODE/transcripts.gtf


# rm *tmp
rm S1_HISAT2_SAM_BAM_FILES/*.sam

mkdir -p REPORT
mv *.summary.txt *.met.txt REPORT

# 5) Predict orfs

TransDecoder.LongOrfs -t transcripts.fa
TransDecoder.Predict -t transcripts.fa --cpu $CPU

exit
