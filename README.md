# Influences of plant genotype, chemotype and environment on the leaf bacterial community

**A. Malacrinò, R. Jakobs, S. Xu, C. Müller**

## Abstract

# Disclaimer

This repository contains the main components used to process the raw data and to analyze it. Raw data is available at NCBI SRA under the BioProject number `XXXXXXXXX`.

Our pipeline included:
* nf-core/ampliseq v2.1.7 [https://github.com/nf-core/ampliseq/](https://github.com/nf-core/ampliseq/)
* R v4.3.2 [https://www.R-project.org/](https://www.R-project.org/)

# Data processing

```bash
nextflow run nf-core/ampliseq -r 2.7.1 -profile singularity \
--input $INDIR \
--FW_primer "AACMGGATTAGATACCCKG" \
--RV_primer "ACGTCATCCCCACCTTCC" \
--outdir $OUTDIR \
--extension "/*_{1,2}.fastq.gz" \
--trunclenf 0 \
--trunclenr 0 \
--skip_qiime \
--skip_barplot \
--skip_abundance_tables \
--skip_alpha_rarefaction \
--skip_diversity_indices \
--skip_ancom \
--max_cpus 16 \
--max_memory '128.GB'

mafft --thread $NTHREADS ASV_seqs.fasta > asv_aligned.fasta

FastTree -gtr -nt < asv_aligned.fasta > tree.tre
```

## Data analysis

Load libraries

```r
library("tidyverse")
library("phyloseq")
library("DESeq2")
library("ggrepel")
library("emmeans")
library("car")
library("lme4")
library("microbiome")
library("ggvenn")
library("decontam")
library("Wrench")
library("RVAideMemoire")
library("picante")
library("decontam")
library("reshape2")
library("ggvenn")
library("MOFA2")
library("RColorBrewer")
library("data.table")
library("psych")
```

```r

```

```r

```

```r

```

```r

```

```r

```

```r

```

```r

```

```r

```

```r

```

```r

```

```r

```

```r

```

```r

```

```r

```

```r

```

```r

```

