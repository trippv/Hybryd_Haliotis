
### Report
```bash
mkdir Report
cp ./*/*.summary.txt Report/
multiqc Report/*
```

### Test if contaminant
```bash
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz

tar -xvf grch38_genome.tar.gz

for i in $(ls *R1_p.fq.gz); do sbatch hisat_align.sh grch38/genome $i; done



# screen
#wget https://github.com/StevenWingett/FastQ-Screen/archive/refs/heads/master.zip
#unzip master.zip 


fastq_screen --aligner bowtie2 --conf fastqscreen.conf --subset 0 --filter 3---- --tag --force --outdir ./fastq_screen/ 012-X04-P68-COI-ICT_S46_L001_R1_001.fastq.fastq
  
# Consider mapping to three genomes (A, B and C), the string '003' produces a file in which reads do not map to genomes A or B, but map (once or more) to genome C.  The string '--1' would generate a file in which reads uniquely map to genome C. Whether reads  map to genome A or B would be ignored.


for i in *qtrim.gz; do fastq_screen --aligner bowtie2 --conf fastqscreen.conf --subset 0 --tag --filter 0000330 --outdir ./fastq_screen/ $i; done
  
# get back the hit_no_genome reads
fastq_screen --conf ../fastqscreen.conf --nohits 012-X04-P68-COI-ICT_S46_L001_R1_001.fastq.headers.tagged.fastq
no_hits=012-X04-P68-COI-ICT_S46_L001_R1_001.fastq.headers.tagged.tagged_filter.fastq


```

