
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

```

# Prediction of protein-coding genes
```bash
See Figure 1.Workflow from this method (Linde et al., 2015; https://doi.org/10.1093/nar/gku1357)
```

# differentially expressed isoforms and Diff splicing events
For the detection of significantly differentially expressed isoforms (DEIs), the mapped reads were utilized with the tools MATS (56), SplicingCompass (57) and DiffSplice (58). Results with a false discovery rate value < 0.05 were considered significant.  (Linde et al., 2015; https://doi.org/10.1093/nar/gku1357)