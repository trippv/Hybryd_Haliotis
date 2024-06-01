```bash
EXPORT=/LUSTRE/apps/bioinformatica/stringtie/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/hisat2/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/diamond_v2.1.8/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3
export PATH=$PATH:$EXPORT


```

Slurm

`hisat2_build.sh`

```bash
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

```

`hisat_align.sh`
```bash
#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH -N 1
#SBATCH --mem=40GB
#SBATCH --ntasks-per-node=20

idx_file=$1
FirstPair=$2 # CFF1_R1_p.fq.gz

CPU=$SLURM_NPROCS

mkdir -p HISAT2_SAM_BAM_FILES

base="${FirstPair%*_*R1_p.fq.gz}"

hisat2  --phred33 --dta -p $CPU \
        -x $idx_file -1 ${base}_R1_p.fq.gz -2 ${base}_R2_p.fq.gz \
        --rg-id=${base} --rg SM:${base} -S HISAT2_SAM_BAM_FILES/${base}.sam \
        --summary-file ${base}.summary.txt --met-file ${base}.met.txt

samtools sort -@ $CPU -o HISAT2_SAM_BAM_FILES/${base}.sorted.bam HISAT2_SAM_BAM_FILES/${base}.sam

exit

```

```
# for file in $(ls *R1_p.fq.gz)
for file in $FirstPair
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_*R1_p.fq.gz}"
hisat2 --phred33 --dta -p $CPU -x $idx_file  \
    -1 ${base}_R1_p.fq.gz -2 ${base}_R2_p.fq.gz \
    --rg-id=${base} --rg SM:${base} -S HISAT2_SAM_BAM_FILES/${base}.sam \
     --summary-file ${base}.summary.txt --met-file ${base}.met.txt
done


```