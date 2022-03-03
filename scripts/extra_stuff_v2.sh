#!/bin/bash
#SBATCH -c 1
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=12:00:00


export PATH=~/omiicsTransfer/software/bedtools2/bin/:$PATH


# v2 improves intron definition by using ucsc table browser introns and intron retention is now only searched in bsj spanning reads, not all reads as before


sample=$1
organism=$2

echo $sample
echo $organism
echo

cd $sample


#if [ organism == "human" ]
#then
### Human
#sample=old_hBr_circ
#genomeSize=/home/mtv/software/IGVTools/genomes/hg19.chrom.sizes
#fa=/home/mtv/faststorage/genomes/human_hg19_july_2010/hg19.fa
#exon_original=/home/mtv/faststorage/genomes/human_hg19_july_2010/gencode.v29lift37.annotation.gffread.exon.bed
#gene=/home/mtv/faststorage/genomes/human_hg19_july_2010/gencode.v29lift37.annotation.gffread.bed
#intron_ucsc=/home/mtv/faststorage/genomes/human_hg19_july_2010/hg19_ucsc_Intron_Gencode_V34lift37.bed
fi
#elif [ organism == "mouse" ]
#then
## Mouse
genomeSize=/home/mtv/software/IGVTools/genomes/mm10.chrom.sizes
fa=/home/mtv/omiicsTransfer/genomes/mm10_GRCm38.p6/GRCm38.p6.genome_simple.fa
exon_original=/home/mtv/omiicsTransfer/genomes/mm10_GRCm38.p6/gencode.vM25.annotation.gffread.exon.merge.bed
gene=/home/mtv/omiicsTransfer/genomes/mm10_GRCm38.p6/gencode.vM25.annotation.gffread.bed
intron_ucsc=/home/mtv/faststorage/genomes/mm10_GRCm38.p6/mm10_ucsc_Intron_Gencode_VM23.bed
# fi

## This reformats the annotation file to include genome region in the name, which makes the output better in the end
cat $exon_original | sed 's/;[[:graph:]]*//g' | awk 'OFS="\t"{print $1,$2,$3,$4"_"$1":"$2"-"$3,0,$6}'  > exons.bed
exon=exons.bed


## Remove empty columns in circRNA candicates file:
mv $sample.circRNA_candidates.annotated.txt OLD.$sample.circRNA_candidates.annotated.txt
grep -v [[:space:]]0[[:space:]][+-][[:space:]]\.[[:space:]]\.[[:space:]]\.[[:space:]]\.[[:space:]]\. OLD.$sample.circRNA_candidates.annotated.txt > $sample.circRNA_candidates.annotated.txt
wc -l OLD.$sample.circRNA_candidates.annotated.txt $sample.circRNA_candidates.annotated.txt


## Finding circRNAs with at least 10 reads that have more than 10% intronic read coverage
printf "internal_circRNA_name\tchr\tstart\tend\tBSJ_reads\tmean_read_coverage\tmean_intron_coverage\tintron_coverage\n" > $sample.circRNA_intron_coverage.txt
cat $sample.circRNA_candidates.annotated.txt | grep -v internal_circRNA_name |  awk '$6>9' | awk 'OFS="\t"{ print $1,$2,$3,$4,$6,$12,$16,$16/$12}' | sort -nrk 5,5 | awk '($7/$6)>0.1' >> $sample.circRNA_intron_coverage.txt
echo "Number of circRNAs with more than 10% intronic coverage"
cat $sample.circRNA_candidates.annotated.txt | grep -v internal_circRNA_name |  awk '$6>9' | awk 'OFS="\t"{ print $1,$2,$3,$4,$6,$12,$16,$16/$12}' | sort -nrk 5,5 | awk '($7/$6)>0.1' | wc -l
echo
## Finding circRNAs with at least 10 reads that have more than 1% intronic read coverage
cat $sample.circRNA_candidates.annotated.txt | grep -v internal_circRNA_name |  awk '$6>9' | awk 'OFS="\t"{ print $1,$2,$3,$4,$6,$12,$16,$16/$12}' | sort -nrk 5,5 | awk '($7/$6)>0.01' >> $sample.circRNA_intron_coverage.1pct.txt
echo "Number of circRNAs with more than 1% intronic coverage"
cat $sample.circRNA_candidates.annotated.txt | grep -v internal_circRNA_name |  awk '$6>9' | awk 'OFS="\t"{ print $1,$2,$3,$4,$6,$12,$16,$16/$12}' | sort -nrk 5,5 | awk '($7/$6)>0.01' | wc -l
echo

echo "These could be either intron retention or due to unannotated exon(s)"
echo
echo "Getting read coverage in introns"
## This was changed in v2 so the intron retention is now only searched in bsj spanning reads, not all reads as before
bedtools coverage -bed -split -a introns.uniq.exon_remove.bed -b $sample.scan.circRNA.sort.bam > $sample.intron.coverage
echo "extracting the introns with >50% of nt covered by reads"
cat $sample.intron.coverage | awk '$10>0.5' > $sample.intron.coverage.50pct
# CircRNAs with more than 1% intronic coverage, if harboring a >50% coverage intron this is shown in last column
printf "internal_circRNA_name\tchr\tstart\tend\tBSJ_reads\tmean_read_coverage\tmean_intron_coverage\tintron_coverage\tcoverage_high_intron\n" > $sample.circRNA_intron_coverage.1pct.50pct_intron.txt
cat $sample.circRNA_intron_coverage.1pct.txt | awk 'OFS="\t"{print $2,$3,$4,$1,$5,$6,$7,$8}' | sortBed | bedtools map -F 1.0 -c 10 -o max -a - -b $sample.intron.coverage.50pct | awk 'OFS="\t"{print $4,$1,$2,$3,$5,$6,$7,$8,$9}'  >> $sample.circRNA_intron_coverage.1pct.50pct_intron.txt


#echo "Cryptic novel exon comparison to annotated exons?"
#
#echo "novel exon sequences"


echo "Map read coverage on intron sequences"
echo "Get introns in circRNAs"
echo "This part was changed in v2"
## Make intron file:
# Getting unique introns that are located in circRNA expression regions
cat $intron_ucsc | awk 'OFS="\t"{print $1,$2,$3,"intron",0,$6}' | sortBed | uniq > introns.uniq.bed
cat $sample.circRNA_candidates.annotated.txt | grep -v internal_circRNA_name | awk 'OFS="\t"{print $2,$3,$4,$1,$6,$7}' | sortBed | bedtools map -c 4 -o distinct -a introns.uniq.bed -b - > $sample.introns.uniq.circ.bed
# Remove introns that overlap an exon from Genncode
bedtools subtract -a introns.uniq.bed -b $exon | sortBed | uniq | awk 'OFS="\t"{print $1,$2,$3,$4"_"$1":"$2"-"$3,1,$6}' | sortBed > introns.uniq.exon_remove.bed
bedtools coverage -a introns.uniq.exon_remove.bed -b $sample.psl.bed > $sample.introns.uniq.exon_remove.coverage.bed
cat $sample.circRNA_candidates.annotated.txt | grep -v internal_circRNA_name | awk 'OFS="\t"{print $2,$3,$4,$1,$6,$7}' | sortBed | bedtools map -c 4 -o distinct -a $sample.introns.uniq.exon_remove.coverage.bed -b - | awk 'OFS="\t"{print $0}' > $sample.introns.uniq.exon_remove.coverage.circ.bed
cat $sample.introns.uniq.exon_remove.coverage.circ.bed | grep circ_ >  $sample.introns.uniq.exon_remove.coverage.onlyCirc.bed
# Mapping novel exons on introns
mapBed -s -F 1.0 -c 4 -o distinct_only -a $sample.introns.uniq.exon_remove.coverage.onlyCirc.bed -b $sample.novel.exons.2reads.bed > $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.bed

# List of all unique introns in circRNA regions:
cat $sample.circRNA_candidates.annotated.txt | grep -v internal_circRNA_name | awk 'OFS="\t"{print $2,$3,$4,$1,$6,$7}' | intersectBed -s -u -a introns.uniq.bed -b - > $sample.all_circRNA_introns.bed
bedtools map -F 1.0 -c 10 -o max -a $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.bed -b $sample.intron.coverage.50pct > $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.intronCov.bed
wc -l introns.bed introns.uniq.bed introns.uniq.exon_remove.bed $sample.introns.uniq.exon_remove.coverage.bed $sample.introns.uniq.exon_remove.coverage.circ.bed $sample.introns.uniq.exon_remove.coverage.onlyCirc.bed $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.bed $sample.all_circRNA_introns.bed

echo
echo "Name of columns in *novelExonMap.bed file"
echo "1: chr, 2: start, 3: end, 4: intron name, 5: filler, 6: strand, 7-10: coverageBed output combared to total reads, 11: contained in circRNA(s), 12: novel exon(s) in intron?"
echo "coverageBed generates theses columns:"
echo "           1) The number of features in B that overlapped the A interval."
echo "           2) The number of bases in A that had non-zero coverage."
echo "           3) The length of the entry in A."
echo "           4) The fraction of bases in A that had non-zero coverage."
echo

#rm introns.bed


# Getting the sequences of novel exons and their phases
bedtools getfasta -fi $fa -fo $sample.novel.exons.2reads.tab -bed $sample.novel.exons.2reads.bed -name -tab -s
cat $sample.novel.exons.2reads.tab | awk '{print $2"XXX"}' | sed 's/.\{3\}/& /g' > temp_p0
cat $sample.novel.exons.2reads.tab | awk '{print $2"XXX"}' | sed 's/^.\(.*\)/\1/' | sed 's/.\{3\}/& /g' > temp_p1
cat $sample.novel.exons.2reads.tab | awk '{print $2"XXX"}' | sed 's/^.\(.*\)/\1/' | sed 's/^.\(.*\)/\1/' | sed 's/.\{3\}/& /g' > temp_p2
paste $sample.novel.exons.2reads.tab temp_p0 temp_p1 temp_p2 | sed 's/(+)//g' | sed 's/(-)//g'  > $sample.novel.exons.2reads.phases.tab




#echo "Liftover"
#### Mouse data
##chain=/home/mtv/ucsc/mm10ToHg19.over.chain
##lift=mm10ToHg19
##comparison_dataset=
#
### Human data
#chain=/home/mtv/ucsc/hg19ToMm10.over.chain
#lift=hg19ToMm10
#comparison_dataset=/home/mtv/faststorage/Karim/old_mBr_circ/old_mBr_circ.circRNA_candidates.annotated.bed
#
##awk 'OFS="\t"{if(NR>1) print $1,$2,$3}' $sample.no_exon_no_circRNA.annotated.bed > no_exon_no_circRNA.annotated.forLiftOver.bed
##liftOver no_exon_no_circRNA.annotated.forLiftOver.bed $chain no_exon_no_circRNA.annotated.$lift.bed no_exon_no_circRNA.annotated.$lift-unmapped.bed
#
#awk 'OFS="\t"{if(NR>1) print $2,$3,$4,$1,$5,$6}' $sample.circRNA_candidates.annotated.txt > $sample.circRNA_candidates.annotated.forLiftover.bed
#awk 'OFS="\t"{if(NR>1) print $1,$2,$2+20,$4,$5,$6}' $sample.circRNA_candidates.annotated.forLiftover.bed > circRNA_candidates.annotated.forLiftover.start.bed
#awk 'OFS="\t"{if(NR>1) print $1,$3-20,$3,$4,$5,$6}' $sample.circRNA_candidates.annotated.forLiftover.bed > circRNA_candidates.annotated.forLiftover.end.bed
#
#liftOver circRNA_candidates.annotated.forLiftover.start.bed $chain circRNA_candidates.annotated.start.$lift.bed circRNA_candidates.annotated.start.$lift.unmapped.bed
#liftOver circRNA_candidates.annotated.forLiftover.end.bed $chain circRNA_candidates.annotated.end.$lift.bed circRNA_candidates.annotated.end.$lift.unmapped.bed
#
#cat circRNA_candidates.annotated.start.$lift.bed circRNA_candidates.annotated.end.$lift.bed | sort -k 4,4 | perl /home/mtv/backup/my-homemade-scripts/Nanopore/combine_liftOver_ends.pl > $sample.circRNA_candidates.annotated.combined.$lift.bed
#bedtools intersect -wo -f 1.0 -F 1.0 -a $sample.circRNA_candidates.annotated.combined.$lift.bed -b $comparison_dataset > $sample.circRNA_candidates.annotated.combined.$lift.human-matches.bed
#wc -l $sample.circRNA_candidates.annotated.combined.$lift.human-matches.bed $sample.circRNA_candidates.annotated.txt $comparison_dataset
#
#rm circRNA_candidates.annotated.forLiftover.start.bed circRNA_candidates.annotated.forLiftover.end.bed circRNA_candidates.annotated.start.*.bed circRNA_candidates.annotated.end.*.bed


