# editing-pipeline

a C to T RNA editing pipeline. Makes use of the Yeo lab [SAILOR](https://github.com/yeolab/sailor) pipeline which they've exported as a singularity image. This pipeline just wraps that image inside snakemake so we can apply it on large datasets and summarise the results. 

The C to T version of SAILOR has been given to the Raj lab by Brian Yee. Thanks Brian!


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

Unfortunately the singularity image can't deal with symlinks so you have to copy these files to references/


