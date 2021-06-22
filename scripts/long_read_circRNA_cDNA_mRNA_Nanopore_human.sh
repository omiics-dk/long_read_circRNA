#!/bin/bash

# Script for running circRNA detaction for Nanopore data


## Variables that need to be adjusted:
organism=human   # Can be human or mouse
sample_name=20191230_20200102_hBrain_mRNA
path_to_data=/home/mtv/faststorage/Karim/cDNA_mRNA_Nanopore


echo
echo "Running script. Data quality filter, mapping and circRNA quantification"
#sh ~/long_read_circRNA_v2/scripts/blat_nanopore_v5.4.sh $path_to_data $sample_name $organism
sbatch --account mtv_AU_project ~/long_read_circRNA_v2/scripts/blat_nanopore_v5.4.sh $path_to_data $sample_name $organism
