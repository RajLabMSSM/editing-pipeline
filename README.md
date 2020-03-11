# editing-pipeline

a C to T RNA editing pipeline. Makes use of the Yeo lab [SAILOR](https://github.com/yeolab/sailor) pipeline which they've exported as a singularity image. This pipeline just wraps that image inside snakemake so we can apply it on large datasets and summarise the results. 

The C to T version of SAILOR has been given to the Raj lab by Brian Yee. Thanks Brian!


## Dependencies

snakemake

## Input files

- list of BAM files (must be aligned with MD flag)
- List of known SNPs to ignore

