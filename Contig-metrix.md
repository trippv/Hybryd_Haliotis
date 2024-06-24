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