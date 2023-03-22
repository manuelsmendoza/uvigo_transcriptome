
RNA sequencing (RNA-Seq) is a powerful technique useful to understand
different biological questions, from gene expression and regulation to
population genetics and evolution. Usually, the standard protocol begins
with the extraction of RNA and its conversion to cDNA in order to be
sequenced using high-throughput technologies. The most often used
application of RNA-Seq is analyzing differential gene expression (DGE),
i.e., comparing the expression of all the genes transcribed under two
different conditions. DGE analysis involved three mainly three steps:

1.  Map the reads for each sample to the reference genome.
2.  Quantify expressed genes and transcripts.
3.  Statistical analysis of gene expression.

> The Sequence Read Archive (SRA, <https://www.ncbi.nlm.nih.gov/sra/>)
> is a public database that stores high-throughput sequencing data to
> hence the open science phylosophy.

## Toy example: Compare the gene expression between human males and females

**Step 1**: Download the reads to analyse the gene expression. The reads
can be downloaded from SRA using `download_reads.R`. The reads
corresponding to chrX from these samples can be found at:
`/mnt/netapp1/share_miba/rnaseq/reads` and the samples information is
stored in the following file
`/mnt/netapp1/share_miba/rnaseq/samples_info.tsv`.

``` bash
Usage: uvigo_transcriptome/download_reads.R [options]

Options:
    -a ACC, --acc=ACC
        Accession number
    -o OUT, --out=OUT
        Output directory
    -c CPU, --cpu=CPU
        Number of threads to use
    -h, --help
        Show this help message and exit
```

``` bash
for RUN in $(awk '{print $1}' samples_info.tsv)
do
  Rscript download_reads.R -a $RUN -o /path/reads_dir -c THREADS
done
```

| Run accession | Sex    |
|:--------------|:-------|
| ERR204889     | Male   |
| ERR204966     | Male   |
| ERR204995     | Male   |
| ERR188273     | Female |
| ERR188060     | Female |
| ERR204896     | Female |

**Step 2**: Download and build the index of the human chromosome X
(chrX) using `get_reference.R`. This scripts can be used to fetch any
chromosome (–seq chr##) or the whole genome (–seq genome) and build
their index.

``` bash
Usage: uvigo_transcriptome/get_reference.R [options]

Options:
    -s SEQ, --seq=SEQ
        Sequence to download
    -o OUT, --out=OUT
        Output file in FASTA fomat (.fasta)
    -a ANN, --ann=ANN
        Download the annotation of the sequence
    -p CDS, --cds=CDS
        Extract protein coding sequences
    -h, --help
        Show this help message and exit
```

``` bash
Rscript get_reference.R -s chrX -o /path/ref_dir -a TRUE -p TRUE
```

**Step 3**: Quantify the gene expression and analyse the differences
between both groups. Use `expression_quant.R` to perform both tasks.

``` bash
Usage: uvigo_transcriptome/expression_quant.R [options]

Options:
    -r REF, --ref=REF
        Reference sequence (.fasta)
    -a ANN, --ann=ANN
        Reference annotation (.gtf)
    -p SPL, --spl=SPL
        Samples information
    -o OUT, --out=OUT
        Output directory
    -t TEC, --tec=TEC
        Type of reference used: genome (geno) or transcriptome (tran)
    -c CPU, --cpu=CPU
        Number of threads to use
    -h, --help
        Show this help message and exit
```

``` bash
Rscript expression_quant.R -r /path/ref_seq.fasta -a /path/ref_ann.gtf -p /path/sample_info.tsv -o /path/out_dir -t geno -c THREADS 
```
