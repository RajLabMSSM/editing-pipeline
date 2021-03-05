# editing-pipeline

Snakemake pipelines for detecting RNA-editing from RNA-seq data.

A collaboration between the Raj lab and the Breen lab.

## Dependencies

snakemake

## Input files

- list of BAM files (must be aligned with MD flag)

- List of known SNPs to ignore
   
I have the latest version of dbSNP as a BED file here:
 /sc/hydra/projects/ad-omics/data/references//hg38_reference/dbSNP/dbSNP_GCF_000001405.38.bed

- human genome build

This is here:

/sc/hydra/projects/ad-omics/data/references//hg38_reference/hg38.fa


## Pipelines included

### AluEditing Index 

Measures global editing rates of Alu elements at all the possible editing types.

Has specific dependencies - to be covered.

### REDItools + Annovar + SNPeff + ... (under construction)

Measures de novo editing genome-wide, applies some basic filters, annotates with gene and transcript consequences

### SAILOR pipeline (not being further developed)

This is a specific C to T RNA editing pipeline. Makes use of the Yeo lab [SAILOR](https://github.com/yeolab/sailor) pipeline which they've exported as a singularity image. This pipeline just wraps that image inside snakemake so we can apply it on large datasets and summarise the results. 

The C to T version of SAILOR has been given to the Raj lab by Brian Yee.



