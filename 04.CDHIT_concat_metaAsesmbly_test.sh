#!/bin/bash

#SBATCH -p cicese
#SBATCH --job-name=CDHIT
#SBATCH --output=cd-hit_metaAssembly-%j.log
#SBATCH --error=cd-hit_metaAssembly-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00



# script para reducir la reduncancia en el ensamble
# Este es un test para evaluar si es posible concatenar y reduci la redundancia con CDHIT

## Previo a este paso, se tiene que cambiar el nombre de los transcritos generados con StringTie

#awk '/^>/ {print ">FF" substr($0, 2)} !/^>/ {print}' FF_cross_transcript.fa > FF_cross_transcript_renamed.fa
#awk '/^>/ {print ">RR_" substr($0, 2)} !/^>/ {print}' /RR_cross_transcript.fa > RR_cross_transcript_renamed.fa
#cat FF_cross_transcript_renamed.fa RR_cross_transcript_renamed.fa > FF_RR_metaAssembly.fa

# rutas
MAIN=/LUSTRE/bioinformatica_data/genomica_funcional/tripp/2024_ROBERTO/ASSEMBLY
transcriptome=$MAIN/metaAssembly/FF_RR_metaAssembly.fa
OUTDIR=$MAIN/metaAssembly/CDHIT
OUTFILE=$OUTDIR/CDHIT_FF_RR_metAssembly.fa
CDHIT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/cdhit/


# correr cdhit con un 99% de similitud (dafault es 90%
# Define CD-HIT-EST parameters
identity=0.95  # 95% identity threshold
word_length=10  # CD-HIT word length parameter
threads=24  # Number of threads for parallel processing

# Run CD-HIT-EST
$CDHIT/cd-hit-est -i "${transcriptome}" -o "${OUTFILE}" -c "${identity}" -n "${word_length}" -T "${threads}" -M 20000

echo "CD-HIT-EST completed. Output saved to ${OUTFILE}")