#!/bin/bash
#SBATCH -c 1
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=12:00:00

source /com/extra/ucsc/2015-04-21/load.sh
# bedtools v2.27.1 fails due to a bug in BEDtools merge (not reporting strand after merge)
source /com/extra/bedtools/2.25.0/load.sh


# v6.1: bug fix and vastly improved method for generation of circRNA_exon_usage.txt

## Human
sample=old_hBr_circ
genomeSize=/home/mtv/software/IGVTools/genomes/hg19.chrom.sizes
fa=/home/mtv/faststorage/genomes/human_hg19_july_2010/hg19.fa
exon_refseq=/home/mtv/faststorage/backup/external_datasets/refFlat_bed/Human_refFlat_exon_hg19_Oct2018.sort.bed
#exon=/home/mtv/faststorage/genomes/human_hg19_july_2010/gencode.v29lift37.annotation.gffread.exon.merge.bed
exon_original=/home/mtv/faststorage/genomes/human_hg19_july_2010/gencode.v29lift37.annotation.gffread.exon.merge.bed
exon_full=/home/mtv/faststorage/genomes/human_hg19_july_2010/gencode.v29lift37.annotation.gffread.exon.bed
circBase=/home/mtv/faststorage/circBase_June-2017/circbase_ucsc-browser/hsa_circRNA_complete.hg19.unique.sort.length.bed
circAtlas=/home/mtv/faststorage/circAtlas/circAtlas2.0_June2019_human_hg19_circRNA.0-based.bed
CIRCpedia=/home/mtv/faststorage/CIRCpedia/CIRCpedia_v2_June2019_human_hg19_All_circRNA.unique.bed

## Mouse
#sample=old_mBr_circ
#genomeSize=/home/mtv/software/IGVTools/genomes/mm10.chrom.sizes
#fa=/home/mtv/faststorage/genomes/mm10/Annotation/Genes/Mus_musculus.GRCm38.87.chr-fix.fa
#exon_refseq=/home/mtv/faststorage/external_datasets/refFlat_bed/Mouse_refFlat_exon_mm10_Oct2018.sort.bed
##exon=/home/mtv/faststorage/genomes/mm10_GRCm38.p6/gencode.vM19.annotation.GeneName-exon.merge.bed
#exon_original=/home/mtv/faststorage/genomes/mm10_GRCm38.p6/gencode.vM19.annotation.exon.bed
#exon_full=/home/mtv/faststorage/genomes/mm10_GRCm38.p6/gencode.vM19.annotation.exon.bed
#circBase=/home/mtv/faststorage/circBase_June-2017/mmu_circRNA_complete.mm10lift.sort.length.bed

## This reformats the annotation file to include genome region in the name, which makes the output better in the end
cat $exon_original | sed 's/;[[:graph:]]*//g' | awk 'OFS="\t"{print $1,$2,$3,$4"_"$1":"$2"-"$3,0,$6}'  > exon_annotation.reformat.bed
exon=exon_annotation.reformat.bed
cat $exon_full | sed 's/;[[:graph:]]*//g' | awk 'OFS="\t"{print $1,$2,$3,$4"_"$1":"$2"-"$3,0,$6}' | sed 's/:/\t/g' | awk 'OFS="\t"{print $1,$2,$3,"exon",$5,$7}' | sort -k 5,5 | sort -k 1 | uniq > exon_specific.reformat.bed
exon_all=exon_specific.reformat.bed


bedtools bed12tobed6 -i $sample.scan.circRNA.psl.bed | awk 'OFS="\t"{print $4,$2,$3,$1,$5,$6}' | grep -v chrM | sortBed > $sample.scan.circRNA.psl.split.bed
bedtools merge -s -d 10 -c 4 -o distinct -i $sample.scan.circRNA.psl.split.bed | awk 'OFS="\t"{print $5,$2,$3,$3-$2,$4,$1}' | sed 's/~/\t/g' | awk 'OFS="\t"{print $1,$2,$3,$6"~"$1"~"$2"~"$3"~"$5,$4,$5}' > $sample.scan.circRNA.psl.split.merge.bed
bedtools flank -g $genomeSize -b 2 -i $sample.scan.circRNA.psl.split.merge.bed > $sample.scan.circRNA.psl.split.merge.flank2.bed
bedtools getfasta -name -fullHeader -fi $fa -bed $sample.scan.circRNA.psl.split.merge.flank2.bed -fo $sample.scan.circRNA.psl.split.merge.flank2.fa
cat  $sample.scan.circRNA.psl.split.merge.flank2.fa | sed ':a;N;$!ba;s/+\n/+\t/g' | sed ':a;N;$!ba;s/-\n/-\t/g' | sed 's/^>//g' | perl /home/mtv/backup/my-homemade-scripts/Nanopore/flank2_combine.pl > $sample.scan.circRNA.psl.split.merge.flank2.fa.bed
cat  $sample.scan.circRNA.psl.split.merge.flank2.fa | sed ':a;N;$!ba;s/+\n/+\t/g' | sed ':a;N;$!ba;s/-\n/-\t/g' | sed 's/^>//g' | perl /home/mtv/backup/my-homemade-scripts/Nanopore/flank2_combine.pl | awk 'OFS="\t"{print $6,$7}' | sort | uniq -c | sort -nrk 1,1 > $sample.scan.circRNA.psl.split.merge.flank2.fa.bed.count
cat $sample.scan.circRNA.psl.split.merge.flank2.fa.bed | grep -P "AGGT$" > $sample.scan.circRNA.psl.split.merge.flank2.AGGT.bed
cat $sample.scan.circRNA.psl.split.merge.flank2.fa.bed | grep -P "ACCT$" > $sample.scan.circRNA.psl.split.merge.flank2.ACCT.bed
cat $sample.scan.circRNA.psl.split.merge.flank2.AGGT.bed | awk 'OFS="\t"{print $1,$2,$3}' | sort -nk 3,3 | sort -nk 2,2 | sort -k 1,1 | uniq -c | awk 'OFS="\t"{print $2,$3,$4,"novelExon",$1,"+"}' | sortBed > $sample.scan.circRNA.psl.split.merge.flank2.posExons.bed
cat $sample.scan.circRNA.psl.split.merge.flank2.ACCT.bed | awk 'OFS="\t"{print $1,$2,$3}' | sort -nk 3,3 | sort -nk 2,2 | sort -k 1,1 | uniq -c | awk 'OFS="\t"{print $2,$3,$4,"novelExon",$1,"-"}' | sortBed > $sample.scan.circRNA.psl.split.merge.flank2.negExons.bed
cat $sample.scan.circRNA.psl.split.merge.flank2.posExons.bed $sample.scan.circRNA.psl.split.merge.flank2.negExons.bed | sortBed | awk 'OFS="\t"{print $1,$2,$3,$5"read_"$4"_"$1":"$2"-"$3,$5,$6}' > $sample.scan.circRNA.psl.split.merge.flank2.allExons.bed
echo
echo "Number of exons"
wc -l $sample.scan.circRNA.psl.split.merge.flank2.posExons.bed $sample.scan.circRNA.psl.split.merge.flank2.negExons.bed $sample.scan.circRNA.psl.split.merge.flank2.allExons.bed
echo
echo "Exons supported by at least 1, 2, 3, 5, 10, 20 and 50 read(s)"
cat $sample.scan.circRNA.psl.split.merge.flank2.allExons.bed | awk 'OFS="\t"{if($5>49) print $0}' > $sample.scan.circRNA.psl.split.merge.flank2.allExons.50reads.bed
cat $sample.scan.circRNA.psl.split.merge.flank2.allExons.bed | awk 'OFS="\t"{if($5>19) print $0}' > $sample.scan.circRNA.psl.split.merge.flank2.allExons.20reads.bed
cat $sample.scan.circRNA.psl.split.merge.flank2.allExons.bed | awk 'OFS="\t"{if($5>9) print $0}' > $sample.scan.circRNA.psl.split.merge.flank2.allExons.10reads.bed
cat $sample.scan.circRNA.psl.split.merge.flank2.allExons.bed | awk 'OFS="\t"{if($5>4) print $0}' > $sample.scan.circRNA.psl.split.merge.flank2.allExons.5reads.bed
cat $sample.scan.circRNA.psl.split.merge.flank2.allExons.bed | awk 'OFS="\t"{if($5>2) print $0}' > $sample.scan.circRNA.psl.split.merge.flank2.allExons.3reads.bed
cat $sample.scan.circRNA.psl.split.merge.flank2.allExons.bed | awk 'OFS="\t"{if($5>1) print $0}' > $sample.scan.circRNA.psl.split.merge.flank2.allExons.2reads.bed
wc -l $sample.scan.circRNA.psl.split.merge.flank2.allExons.bed $sample.scan.circRNA.psl.split.merge.flank2.allExons.2reads.bed $sample.scan.circRNA.psl.split.merge.flank2.allExons.3reads.bed $sample.scan.circRNA.psl.split.merge.flank2.allExons.5reads.bed $sample.scan.circRNA.psl.split.merge.flank2.allExons.10reads.bed $sample.scan.circRNA.psl.split.merge.flank2.allExons.20reads.bed $sample.scan.circRNA.psl.split.merge.flank2.allExons.50reads.bed

echo
echo "Number of exons not in Gencode (with 95% similarity)"
bedtools intersect -v -f 0.95 -F 0.95 -a $sample.scan.circRNA.psl.split.merge.flank2.allExons.bed -b $exon > $sample.scan.circRNA.psl.split.merge.flank2.allExons.notGencode.bed
bedtools intersect -v -f 0.95 -F 0.95 -a $sample.scan.circRNA.psl.split.merge.flank2.allExons.2reads.bed -b $exon > $sample.novel.exons.2reads.bed
bedtools intersect -v -f 0.95 -F 0.95 -a $sample.novel.exons.2reads.bed -b $exon_refseq > $sample.novel.exons.2reads.filter00.bed
bedtools intersect -v -f 0.95 -F 0.95 -a $sample.novel.exons.2reads.filter00.bed -b $exon_all > $sample.novel.exons.2reads.filter0.bed

mv $sample.novel.exons.2reads.filter0.bed $sample.novel.exons.2reads.filter.bed


wc -l $sample.scan.circRNA.psl.split.merge.flank2.allExons.notGencode.bed $sample.novel.exons.2reads.bed $sample.novel.exons.2reads.filter00.bed $sample.novel.exons.2reads.filter0.bed $sample.novel.exons.2reads.filter.bed
rm $sample.novel.exons.2reads.filter0.bed $sample.novel.exons.2reads.filter00.bed

echo
echo "Getting 2read circRNA info for use with NMD test"
intersectBed -wo -a $sample.novel.exons.2reads.filter.bed -b $exon | awk 'OFS="\t"{print $1,$2,$3,$4,$7,$8,$9,$10,$3-$2,$9-$8,($9-$8)-($3-$2),(($9-$8)-($3-$2))/3}' > $sample.novel.cryptic.spliced.exons.txt


# Combine novel exons with atleast 2 supporting reads (not found in Gencode) with all Gencode exons
cat $exon $sample.novel.exons.2reads.filter.bed | sortBed | uniq | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6}' > New_exon_and_Gencode_exon.bed

bedtools map -split -c 4 -o distinct -a $sample.scan.circRNA.psl.bed -b New_exon_and_Gencode_exon.bed > $sample.scan.circRNA.psl.circRNA-exons.bed
bedtools map -c 4 -o distinct -a $sample.scan.circRNA.psl.circRNA-exons.bed -b New_exon_and_Gencode_exon.bed > $sample.scan.circRNA.psl.genomic-exons.bed

# Prepare for comparison
cat $sample.scan.circRNA.psl.genomic-exons.bed | awk 'OFS="\t"{print $1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$4}' | sed 's/~/\t/g'| awk 'OFS="\t"{print $14,$2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | sortBed > temp.genomic-exons.bed

cat $sample.scan.circRNA.psl.annot.combine.txt | sortBed > $sample.scan.circRNA.psl.annot.combine.sort.txt

cat $sample.circRNA_candidates.annotated.txt | awk 'NR>1,OFS="\t"{print $2,$3,$4,$1,$6,$7,$10}' | sortBed > circRNA_candidates.annotated.txt.forExonAnalysis.bed
bedtools map -f 0.95 -F 0.95 -c 4,7 -o distinct -a $sample.scan.circRNA.psl.annot.combine.sort.txt -b circRNA_candidates.annotated.txt.forExonAnalysis.bed | awk 'OFS="\t"{print $4,$2,$3,$1,0,"+",$10}' | sed 's/,circ/\t/g' | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7}' | sortBed > reads.annot.temp.bed

bedtools map -c 7 -o distinct -a temp.genomic-exons.bed -b reads.annot.temp.bed | awk 'OFS="\t"{print $4,$2,$3,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' | sortBed > $sample.scan.circRNA.psl.genomic-exons.annot.bed

rm reads.annot.temp.bed temp.genomic-exons.bed $sample.scan.circRNA.psl.split.merge.flank2.posExons.bed $sample.scan.circRNA.psl.split.merge.flank2.negExons.bed
rm $sample.scan.circRNA.psl.annot.combine.sort.txt




echo "Starting the list with number of used exons for each circRNA"
date
mkdir /home/mtv/temp/job.$SLURM_JOBID
cat $sample.scan.circRNA.psl.genomic-exons.annot.bed | awk 'OFS="\t"{print $13,$14,$15}' | sort -k 1,1 -T /home/mtv/temp/job.$SLURM_JOBID > temp1
cat temp1 | sort -k 2,2 -T /home/mtv/temp/job.$SLURM_JOBID > temp2
cat temp2 | sort -k 3,3 -T /home/mtv/temp/job.$SLURM_JOBID | uniq -c > $sample.scan.circRNA.psl.genomic-exons.annot.uniq.bed
rm -r temp1 temp2 /home/mtv/temp/job.$SLURM_JOBID




## New major improvement in version 6.1:
input=$sample.scan.circRNA.psl.genomic-exons.annot.uniq.bed

touch circRNA_exon_usage.txt
rm circRNA_exon_usage.txt

#cat $input | awk '{print $4}' | sort | uniq | awk '{if ($1 != ".") print $1}'> circRNA-list

while IFS='' read -r circRNA || [[ -n "$circRNA" ]]; do
        ## Getting the exons the specific circRNA checked now. Both novel and annotated exons are included, but  novel only if they have 50 or more reads.

       ### If running on a system with SLURM queueing system (This will take a long time):
#sbatch core_circ_grep_10reads.sh $circRNA $input $sample

	### If no SLURM queueing system:
        #Note: This removes novel exons with fewer than 10 reads.
        grep $circRNA $input | awk '{print $3}' | sed 's/,/\n/g' | sort | uniq | grep -v "^[123456789]read_novelExon" > temp_exon_list
        echo
        echo "circRNA: "$circRNA
        echo "exons in circRNA"
        cat temp_exon_list
        echo

        while IFS='' read -r exon || [[ -n "$exon" ]]; do
                exon_hit=$(grep $circRNA $SLURM_SUBMIT_DIR/$sample.scan.circRNA.psl.genomic-exons.annot.uniq.bed | awk '{print $1,$2}' | grep $exon | awk '{split($0,a,","); sum += a[1]} END {print sum}')
                circRNA_coverage=$(grep $circRNA $SLURM_SUBMIT_DIR/$sample.scan.circRNA.psl.genomic-exons.annot.uniq.bed | awk '{print $1,$3}' | grep $exon | awk '{split($0,a,","); sum += a[1]} END {print sum}')
                printf "$circRNA\t$exon_hit\t$circRNA_coverage\t$exon" | awk 'OFS="\t"{if($3 != 0) {print $1,$2,$3,$2/$3,$4;}}' >> $SLURM_SUBMIT_DIR/$circRNA.circRNA_exon_usage.txt

        done < temp_exon_list



done < circRNA-list



#cat circRNA_exon_usage.txt | awk '$2>9' | awk '$4<0.9' | awk '$4>0.1' > circRNA_alternative_exon_usage.txt
#wc -l circRNA_exon_usage.txt circRNA_alternative_exon_usage.txt
#
#mkdir DELETE
#mv core_circ_grep.sh-*.out DELETE
#mv *.circRNA_exon_usage.txt DELETE

echo Done
date
