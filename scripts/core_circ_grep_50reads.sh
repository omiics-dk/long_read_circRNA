#!/bin/bash
#SBATCH -c 1
#SBATCH -p normal
#SBATCH --mem=12g
#SBATCH --time=12:00:00

source /com/extra/ucsc/2015-04-21/load.sh
# bedtools v2.27.1 fails due to a bug in BEDtools merge (not reporting strand after merge)
source /com/extra/bedtools/2.25.0/load.sh


# change directory to the local scratch-directory, and run:
cd /scratch/$SLURM_JOBID

circRNA=$1
input=$2
sample=$3

######################### This is core_circ_grep_50reads'.sh
	#Note: This removes novel exons with fewer than 50 reads. Only for panel datasets. Otherwise use cut off 10
        grep $circRNA $SLURM_SUBMIT_DIR/$input | awk '{print $3}' | sed 's/,/\n/g' | sort | uniq | grep -v "^[123456789]read_novelExon" | grep -v "^[123][1234567890]read_novelExon" | grep -v "^4[123456789]read_novelExon" > temp_exon_list
	echo
	echo "circRNA: "$circRNA
	echo "exons in circRNA"
	cat temp_exon_list
	echo


        while IFS='' read -r exon || [[ -n "$exon" ]]; do
                exon_hit=$(grep $circRNA $SLURM_SUBMIT_DIR/$sample.scan.circRNA.psl.genomic-exons.annot.uniq.bed | awk '{print $1,$2}' | grep $exon | awk '{split($0,a,","); sum += a[1]} END {print sum}')
                circRNA_coverage=$(grep $circRNA $SLURM_SUBMIT_DIR/$sample.scan.circRNA.psl.genomic-exons.annot.uniq.bed | awk '{print $1,$3}' | grep $exon | awk '{split($0,a,","); sum += a[1]} END {print sum}')
#                printf "$circRNA\t$exon_hit\t$circRNA_coverage\t$exon" | awk 'OFS="\t"{print $1,$2,$3,$2/$3,$4}' > $SLURM_SUBMIT_DIR/$circRNA.circRNA_exon_usage.txt
		printf "$circRNA\t$exon_hit\t$circRNA_coverage\t$exon" | awk 'OFS="\t"{if($3 != 0) {print $1,$2,$3,$2/$3,$4;}}' >> $SLURM_SUBMIT_DIR/$circRNA.circRNA_exon_usage.txt

        done < temp_exon_list
######################### core_circ_grep.sh done

