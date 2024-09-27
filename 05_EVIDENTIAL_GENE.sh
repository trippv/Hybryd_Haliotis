#!/bin/bash

#SBATCH -p cicese
#SBATCH --job-name=EviGene
#SBATCH --output=EvidentialGene_traacds-%j.log
#SBATCH --error=EvidentialGene_traacds-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

# Descripcion: Script para reducir el numero de transcritos de un ensamble 
# about: http://arthropods.eugenes.org/EvidentialGene/evigene/docs/EvidentialGene_howto.txt


# cargar modulos
#export some soft
export PATH=$PATH:/LUSTRE/apps/bioinformatica/bowtie2/bin/
export PATH=$PATH:/LUSTRE/apps/bioinformatica/samtools-1.17/bin/
export PATH=$PATH:/LUSTRE/apps/bioinformatica/jellyfish-2.3.0/
export PATH=$PATH:/LUSTRE/apps/bioinformatica/salmon/bin/

export PATH=$PATH:/LUSTRE/apps/bioinformatica/ncbi-blast-2.13.0+/bin/blastn
export PATH=$PATH:/LUSTRE/apps/bioinformatica/ncbi-blast-2.13.0+/bin/makeblastdb
export PATH=$PATH:/LUSTRE/bioinformatica_data/genomica_funcional/apps/exonerate-2.2.0/src/util/fastanrdb
export PATH=$PATH:/LUSTRE/bioinformatica_data/genomica_funcional/apps/cdhit/cd-hit-est
export PATH=$PATH:/LUSTRE/bioinformatica_data/genomica_funcional/apps/cdhit/cd-hit

module load spack
module load blast-plus/2.13.0-gcc-12.2.0-k4ltzsw 

# evigene path
EVIGENE=/LUSTRE/bioinformatica_data/genomica_funcional/apps/evigene/scripts/prot/tr2aacds.pl
FORMAT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/evigene/scripts/rnaseq/trformat.pl

# rutas
MAIN=/LUSTRE/bioinformatica_data/genomica_funcional/tripp/2024_ROBERTO
EVIGENE_DIR=$MAIN/EVIGENE

# merged assembly
transcriptome=$EVIGENE_DIR/MetaAssembly.fa # ensamble concatenado de RR y FF (Hisat-stringtie) y RF (RNAspades)

# crear el folder con el indice en caso de que no exista
if [ ! -d "$EVIGENE_DIR" ]; then
    mkdir "$EVIGENE_DIR"
    echo "Folder '$EVIGENE_DIR' created."
else
    echo "Folder '$EVIGENE_DIR' already exists."
fi


cd $EVIGENE_DIR

### run evigene
$EVIGENE -cdnaseq $transcriptome -NCPU=24 -MAXMEM=100000 -MINAA=100



