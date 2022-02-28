#!/bin/bash

# Script for running circRNA detaction for Nanopore data


## Variables that need to be adjusted:
organism="mouse"   # Can be human or mouse


for sample in a #b c d e f
do
if [ $sample == "a" ]
then
sample_name=E14_1
path_to_data=~/primary_data/backup/Karim/Matejas_Nanopore_seq/ONT_circRNA/fastq/pass/barcode01/

elif [ $sample == "b" ]
then
sample_name=E14_3
path_to_data=~/primary_data/backup/Karim/Matejas_Nanopore_seq/ONT_circRNA/fastq/pass/barcode02/

elif [ $sample == "c" ]
then
sample_name=E14_4
path_to_data=~/primary_data/backup/Karim/Matejas_Nanopore_seq/ONT_circRNA/fastq/pass/barcode03/

elif [ $sample == "d" ]
then
sample_name=P0_5
path_to_data=~/primary_data/backup/Karim/Matejas_Nanopore_seq/ONT_circRNA/fastq/pass/barcode04/

elif [ $sample == "e" ]
then
sample_name=P0_7
path_to_data=~/primary_data/backup/Karim/Matejas_Nanopore_seq/ONT_circRNA/fastq/pass/barcode05/

elif [ $sample == "f" ]
then
sample_name=P0_9
path_to_data=~/primary_data/backup/Karim/Matejas_Nanopore_seq/ONT_circRNA/fastq/pass/barcode06/
fi

echo
echo "Running script. Data quality filter, mapping and circRNA quantification"
sbatch --account mtv_AU_project ~/omiicsTransfer/long_read_circRNA_v2/scripts/blat_nanopore_v5.5.sh $path_to_data $sample_name $organism
sbatch --account mtv_AU_project ~/omiicsTransfer/long_read_circRNA_v2/scripts/novel_exons_and_alternative_usage_v7.0.sh $sample_name $organism
#sbatch --account mtv_AU_project temp_extra.sh $sample_name $organism
done
