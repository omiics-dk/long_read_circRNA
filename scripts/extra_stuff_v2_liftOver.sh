#!/bin/bash

source /com/extra/ucsc/2015-04-21/load.sh


# v2 improves intron definition by using ucsc table browser introns and intron retention is now only searched in bsj spanning reads, not all reads as before


## Human
sample=old_hBr_circ
genomeSize=/home/mtv/software/IGVTools/genomes/hg19.chrom.sizes
fa=/home/mtv/faststorage/genomes/human_hg19_july_2010/hg19.fa
exon_original=/home/mtv/faststorage/genomes/human_hg19_july_2010/gencode.v29lift37.annotation.gffread.exon.bed
gene=/home/mtv/faststorage/genomes/human_hg19_july_2010/gencode.v29lift37.annotation.gffread.bed
circBase=/home/mtv/faststorage/circBase_June-2017/circbase_ucsc-browser/hsa_circRNA_complete.hg19.unique.sort.length.bed
circAtlas=/home/mtv/faststorage/circAtlas/circAtlas2.0_June2019_human_hg19_circRNA.0-based.bed
CIRCpedia=/home/mtv/faststorage/CIRCpedia/CIRCpedia_v2_June2019_human_hg19_All_circRNA.unique.bed
intron_ucsc=/home/mtv/faststorage/genomes/human_hg19_july_2010/hg19_ucsc_Intron_Gencode_V34lift37.bed

## Mouse
#sample=old_mBr_circ
#genomeSize=/home/mtv/software/IGVTools/genomes/mm10.chrom.sizes
#fa=/home/mtv/faststorage/genomes/mm10/Annotation/Genes/Mus_musculus.GRCm38.87.chr-fix.fa
#exon_original=/home/mtv/faststorage/genomes/mm10_GRCm38.p6/gencode.vM19.annotation.exon.bed
#gene=/home/mtv/faststorage/genomes/mm10_GRCm38.p6/gencode.vM19.annotation.bed
#circBase=/home/mtv/faststorage/circBase_June-2017/mmu_circRNA_complete.mm10lift.sort.length.bed
#intron_ucsc=/home/mtv/faststorage/genomes/mm10_GRCm38.p6/mm10_ucsc_Intron_Gencode_VM23.bed

## This reformats the annotation file to include genome region in the name, which makes the output better in the end
cat $exon_original | sed 's/;[[:graph:]]*//g' | awk 'OFS="\t"{print $1,$2,$3,$4"_"$1":"$2"-"$3,0,$6}'  > exons.bed
exon=exons.bed





echo "Liftover"
### Mouse data
#chain=/home/mtv/ucsc/mm10ToHg19.over.chain
#lift=mm10ToHg19
#comparison_dataset=/home/mtv/faststorage/Karim/old_hBr_circ/old_hBr_circ.circRNA_candidates.annotated.bed

## Human data
chain=/home/mtv/ucsc/hg19ToMm10.over.chain
lift=hg19ToMm10
comparison_dataset=/home/mtv/faststorage/Karim/old_mBr_circ/old_mBr_circ.circRNA_candidates.annotated.bed

#awk 'OFS="\t"{if(NR>1) print $1,$2,$3}' $sample.no_exon_no_circRNA.annotated.bed > no_exon_no_circRNA.annotated.forLiftOver.bed
#liftOver no_exon_no_circRNA.annotated.forLiftOver.bed $chain no_exon_no_circRNA.annotated.$lift.bed no_exon_no_circRNA.annotated.$lift-unmapped.bed

awk 'OFS="\t"{if(NR>1) print $2,$3,$4,$1,$5,$6}' $sample.circRNA_candidates.annotated.txt > $sample.circRNA_candidates.annotated.forLiftover.bed
awk 'OFS="\t"{if(NR>1) print $1,$2,$2+20,$4,$5,$6}' $sample.circRNA_candidates.annotated.forLiftover.bed > circRNA_candidates.annotated.forLiftover.start.bed
awk 'OFS="\t"{if(NR>1) print $1,$3-20,$3,$4,$5,$6}' $sample.circRNA_candidates.annotated.forLiftover.bed > circRNA_candidates.annotated.forLiftover.end.bed

liftOver circRNA_candidates.annotated.forLiftover.start.bed $chain circRNA_candidates.annotated.start.$lift.bed circRNA_candidates.annotated.start.$lift.unmapped.bed
liftOver circRNA_candidates.annotated.forLiftover.end.bed $chain circRNA_candidates.annotated.end.$lift.bed circRNA_candidates.annotated.end.$lift.unmapped.bed

cat circRNA_candidates.annotated.start.$lift.bed circRNA_candidates.annotated.end.$lift.bed | sort -k 4,4 | perl /home/mtv/backup/my-homemade-scripts/Nanopore/combine_liftOver_ends.pl > $sample.circRNA_candidates.annotated.combined.$lift.bed
bedtools intersect -wo -f 1.0 -F 1.0 -a $sample.circRNA_candidates.annotated.combined.$lift.bed -b $comparison_dataset > $sample.circRNA_candidates.annotated.combined.$lift.human-matches.bed
wc -l $sample.circRNA_candidates.annotated.combined.$lift.human-matches.bed $sample.circRNA_candidates.annotated.txt $comparison_dataset

rm circRNA_candidates.annotated.forLiftover.start.bed circRNA_candidates.annotated.forLiftover.end.bed circRNA_candidates.annotated.start.*.bed circRNA_candidates.annotated.end.*.bed


