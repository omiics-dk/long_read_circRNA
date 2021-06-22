#!/bin/bash

# Script for running circRNA detaction for Nanopore data


## Variables that need to be adjusted:
organism=human   # Can be human or mouse
sample_name=human_brain_100k
path_to_data=~/long_read_circRNA/test_fastq

## Set type to:
## "panel" for data sequenced from a limited panel of circRNAs at high depth or
## "nick" for data sequenced from complete circRNA enriched pool with partial nicking
type=nick



echo
echo "Running script. Data quality filter, mapping and circRNA quantification"
sh ~/long_read_circRNA/scripts/blat_nanopore_v5.4.sh $path_to_data $sample_name $organism
### Alternatively, if running on a cluster with the SLURM queueing system you can use the following command:
#sbatch ~/long_read_circRNA/scripts/blat_nanopore_v5.4.sh $path_to_data $sample_name $organism
