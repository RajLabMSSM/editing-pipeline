#!/bin/bash
#BSUB -J "AEI"
#BSUB -P acc_buxbaj01a
#BSUB -q premium
#BSUB -o ./aei_%J.%I
#BSUB -u enrico.mossotto@mssm.edu
#BSUB -L /bin/bash
#BSUB -W 05:00
#BSUB -n 10
#BSUB -R rusage[mem=4000]
#BSUB -R span[hosts=1]
#BSUB -M 20GB

# references stored here:
/sc/arion/projects/breen_lab/AEI
module load rnaeditingindexer

. ./variables

DIRE=$(dirname ${BAM})
BB=$(basename ${BAM})

RNAEditingIndex -d $DIRE -f $BB -o . --genes_expression /sc/hydra/scratch/mossoe02/Resources/GenesExpression/HomoSapiens/ucscHg38GTExGeneExpression.bed.gz --refseq /sc/hydra/scratch/mossoe02/Resources/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz --snps /sc/hydra/scratch/mossoe02/Resources/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz -gf /sc/hydra/scratch/mossoe02/Resources/Genomes/HomoSapiens/ucscHg38Genome.fa -rb /sc/hydra/scratch/mossoe02/Resources/Regions/HomoSapiens/ucscHg38Alu.bed.gz --genome UserProvided  --paired_end --stranded

rm *.cmpileup
