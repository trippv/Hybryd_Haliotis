# Hybryd_Haliotis


# Haliotis rufescens
/LUSTRE/bioinformatica_data/genomica_funcional/Oyervides/data/OA_abulon/RNASEQ/Genomes/Haliotis_rufescens

# Haliotis fulgens
/LUSTRE/bioinformatica_data/genomica_funcional/Oyervides/data/OA_abulon/RNASEQ/Genomes/Haliotis_fulgens/v2


# Red
awk '{print $3}' *.gtf| sort | uniq -c
 591736 CDS
 665954 exon
  77553 five_prime_utr
  42838 gene
     22 Selenocysteine
  55645 start_codon
  55603 stop_codon
  54432 three_prime_utr
  70020 transcript

# Green
  awk '{print $3}' *.gff3 | sort | uniq -c
 324879 CDS
 329100 exon
  42489 gene
  45161 mRNA
     96 ncRNA
     95 ncRNA_gene