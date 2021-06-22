# long_read_circRNA
CircRNA detection in nanopore data

#################################

Info for using the scripts in long_read_circRNA

Scripts are made for linux and long_read_circRNA directory should be placed directly in the home directory.

The scripts supplied here need the following software so work
Version numbers are shown, but other versions of these programs might also work:
bedtools v2.29.2
minimap2 v2.17
nanofilt v2.6.0
 

The directory long_read_circRNA has three subdirectories:
data - Holds all the data files needed for analysis
scripts - Has all the scripts used
test_fastq - Has a test dataset with human nanopore circRNA fastq data which can be used to test the scripts


To run nanopore circRNA detection the following script needs to be updated with the relevant info, such as
sample name, path to data and organism. The script is currently made ready to analyze the provided test data:
~/long_read_circRNA/scripts/long_read_circRNA.sh

To start the circRNA quantification:
sh ~/long_read_circRNA/scripts/long_read_circRNA.sh


The scripts output a number of files, the primary one being:
[sample].circRNA_candidates.annotated.txt

This file shows all circRNAs detected in the nanopore data, listed with the highest expressed at the top.
