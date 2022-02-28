#!/bin/bash

# Script for running circRNA detaction for Nanopore data


#The Mateja's Nanopore data is uploaded now to the cluster. It should be around 180 GB including 2,324 fasta5 files.
#The barcodes are numbers 1-6 from PCR-cDNA barcoding kit (SQK-PCB109).
#
#Barcodes 1-3 are three replicates from E14 time point (samples 1, 3, 4).
#Barcodes 4-6 are three replicates from P0 time point (samples 5, 7, 9).


## Variables that need to be adjusted:
organism=mouse   # Can be human or mouse
sample_name=E14_1
path_to_data=~/primary_data/backup/Karim/Matejas_Nanopore_seq/ONT_circRNA/fastq/pass/barcode01/

## Set type to:
## "panel" for data sequenced from a limited panel of circRNAs at high depth or
## "nick" for data sequenced from complete circRNA enriched pool with partial nicking
type=panel



echo
echo "Running script. Data quality filter, mapping and circRNA quantification"
sh ~/long_read_circRNA_v2/scripts/blat_nanopore_v5.4.sh $path_to_data $sample_name $organism
### Alternatively, if running on a cluster with the SLURM queueing system you can use the following command:
#sbatch ~/long_read_circRNA/scripts/blat_nanopore_v5.4.sh $path_to_data $sample_name $organism
