Transrate easly calculate contig metrics for an assembly
```bash

EXPORT=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin
export PATH=$PATH:$EXPORT

export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH


transrate --assembly denovo-transcript.fa --threads 10
transrate --assembly refbased-transcript.fa --threads 10
transrate --assembly rnaspades.fasta --threads 10
transrate --assembly Trinity.fasta --threads 10

```

Busco

```bash
/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Haliotis/03.Merges/REFBASED_MODE/BUSCO
mkdir summaries
cp ./run_*odb*/short_summary* summaries

module load R-4.3.1

EXPORT=/home/rgomez/bin/busco-master/scripts/
export PATH=$EXPORT:$PATH

generate_plot.py --working_directory ./summaries/

cp summaries/busco_figure.R .
sed -i 's/_odb10//g' summaries/busco_figure.R
Rscript summaries/busco_figure.R
# firefox summaries/busco_figure.png
cp summaries/busco_figure.R .
```

## gffcompare
Gffcompare can be used to evaluate and compare the accuracy of transcript assemblers - in terms of their structural correctness (exon/intron coordinates). 

In order to compare the baseline de novo transcript assembly accuracy, for instance, both different assemblers, should be run without using any reference annotation data (i.e. no -G or -g options were used).


```bash
 EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/gffcompare
export PATH=$PATH:$EXPORT

cd /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Haliotis/03.Merges/REFBASED_MODE

# gffcompare -r hybrid.gtf -i stringtie_gtf_list.txt

gffcompare -r hybrid.gtf transcripts.gtf

cd /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Haliotis/03.Merges/DENOVO_MODE/

gffcompare -r hybrid.gtf transcripts.gtf



cd /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Haliotis/Contig-metrics/rnaspades-hisat-alignment/HISAT2_SAM_BAM_FILES

gffcompare -r hybrid.gtf transcripts.gtf
# resulted in zero cause hybrid.gtf does not match to ids from transcripts.gtf

# The measures of “sensitivity” and “precision” are calculated at various levels (nucleotide, exon, intron, transcript, gene) for each input file and reported in this .stats file.

```