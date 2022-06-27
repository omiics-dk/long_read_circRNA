#!/bin/bash
#SBATCH -c 1
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=12:00:00


sample=$1
organism=$2
reference_path=$3
scriptFolder=$4

echo $sample
echo $organism
echo
date

cd $5

# v6.1: bug fix and vastly improved method for generation of circRNA_exon_usage.txt
# v7.0: Added the script previously called extra_stuff_v2.sh to this script, so now intron analyses, and phasing is done here

if [ $organism == "human" ]
then
### Human
genomeSize=$reference_path/human/hg19.chrom.sizes
fa=$reference_path/human/hg19.fa
exon_refseq=$reference_path/human/Human_refFlat_exon_hg19_Oct2018.sort.bed
exon_original=$reference_path/human/gencode.v37lift37.annotation.gffread.exon.merge.bed
exon_full=$reference_path/human/gencode.v37lift37.annotation.gffread.exon.bed
intron_ucsc=$reference_path/human/hg19_ucsc_Intron_Gencode_V34lift37.bed
elif [ $organism == "mouse" ]
then
### Mouse
genomeSize=$reference_path/mouse/mm10.chrom.sizes
fa=$reference_path/mouse/Mus_musculus.GRCm38.87.chr-fix.fa
exon_refseq=$reference_path/mouse/Mouse_refFlat_exon_mm10_Oct2018.sort.bed
exon_original=$reference_path/mouse/gencode.vM25.annotation.gffread.exon.merge.bed
exon_full=$reference_path/mouse/gencode.vM25.annotation.gffread.exon.bed
intron_ucsc=$reference_path/mouse/mm10_ucsc_Intron_Gencode_VM23.bed
fi

echo "Reformating annotation files"
date

## This reformats the annotation file to include genome region in the name, which makes the output better in the end
cat $exon_original | sed 's/;[[:graph:]]*//g' | awk 'OFS="\t"{print $1,$2,$3,$4"_"$1":"$2"-"$3,0,$6}'  > exon_annotation.reformat.bed
date
exon=exon_annotation.reformat.bed
cat $exon_full | sed 's/;[[:graph:]]*//g' | awk 'OFS="\t"{print $1,$2,$3,$4"_"$1":"$2"-"$3,0,$6}' | sed 's/:/\t/g' | awk 'OFS="\t"{print $1,$2,$3,"exon",$5,$7}' | sort -k 5,5 | sort -k 1 | uniq > exon_specific.reformat.bed
exon_all=exon_specific.reformat.bed

date
echo
echo "Producing summary of exons"

bedtools bed12tobed6 -i $sample.scan.circRNA.psl.bed | awk 'OFS="\t"{print $4,$2,$3,$1,$5,$6}' | grep -v chrM | sortBed > $sample.scan.circRNA.psl.split.bed
bedtools merge -s -d 10 -c 4,6 -o distinct -i $sample.scan.circRNA.psl.split.bed | awk 'OFS="\t"{print $4,$2,$3,$3-$2,$1,$5}' | sed 's/~/\t/g' | awk 'OFS="\t"{print $1,$2,$3,$5"~"$1"~"$2"~"$3"~"$7,$4,$5}' > $sample.scan.circRNA.psl.split.merge.bed
#cat $genomeSize
bedtools flank -g $genomeSize -b 2 -i $sample.scan.circRNA.psl.split.merge.bed > $sample.scan.circRNA.psl.split.merge.flank2.bed
bedtools getfasta -nameOnly -fi $fa -bed $sample.scan.circRNA.psl.split.merge.flank2.bed -fo $sample.scan.circRNA.psl.split.merge.flank2.fa
cat $sample.scan.circRNA.psl.split.merge.flank2.fa | sed ':a;N;$!ba;s/+\n/+\t/g' | sed ':a;N;$!ba;s/-\n/-\t/g' | sed 's/^>//g' | perl $scriptFolder/flank2_combine.pl > $sample.scan.circRNA.psl.split.merge.flank2.fa.bed
cat $sample.scan.circRNA.psl.split.merge.flank2.fa | sed ':a;N;$!ba;s/+\n/+\t/g' | sed ':a;N;$!ba;s/-\n/-\t/g' | sed 's/^>//g' | perl $scriptFolder/flank2_combine.pl | awk 'OFS="\t"{print $6,$7}' | sort | uniq -c | sort -nrk 1,1 > $sample.scan.circRNA.psl.split.merge.flank2.fa.bed.count
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
bedtools intersect -nonamecheck -v -f 0.95 -F 0.95 -a $sample.scan.circRNA.psl.split.merge.flank2.allExons.bed -b $exon > $sample.scan.circRNA.psl.split.merge.flank2.allExons.notGencode.bed
bedtools intersect -nonamecheck -v -f 0.95 -F 0.95 -a $sample.scan.circRNA.psl.split.merge.flank2.allExons.2reads.bed -b $exon > $sample.novel.exons.2reads.bed
bedtools intersect -nonamecheck -v -f 0.95 -F 0.95 -a $sample.novel.exons.2reads.bed -b $exon_refseq > $sample.novel.exons.2reads.filter00.bed
bedtools intersect -nonamecheck -v -f 0.95 -F 0.95 -a $sample.novel.exons.2reads.filter00.bed -b $exon_all > $sample.novel.exons.2reads.filter0.bed

mv $sample.novel.exons.2reads.filter0.bed $sample.novel.exons.2reads.filter.bed


wc -l $sample.scan.circRNA.psl.split.merge.flank2.allExons.notGencode.bed $sample.novel.exons.2reads.bed $sample.novel.exons.2reads.filter00.bed $sample.novel.exons.2reads.filter.bed
rm $sample.novel.exons.2reads.filter00.bed

echo
echo "Getting 2read circRNA info for use with NMD test"
intersectBed -nonamecheck -wo -a $sample.novel.exons.2reads.filter.bed -b $exon | awk 'OFS="\t"{print $1,$2,$3,$4,$7,$8,$9,$10,$3-$2,$9-$8,($9-$8)-($3-$2),(($9-$8)-($3-$2))/3}' > $sample.novel.cryptic.spliced.exons.txt


# Combine novel exons with atleast 2 supporting reads (not found in Gencode) with all Gencode exons
cat $exon $sample.novel.exons.2reads.filter.bed | sortBed | uniq | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6}' > New_exon_and_Gencode_exon.bed

bedtools map -nonamecheck -split -c 4 -o distinct -a $sample.scan.circRNA.psl.bed -b New_exon_and_Gencode_exon.bed > $sample.scan.circRNA.psl.circRNA-exons.bed
bedtools map -nonamecheck -c 4 -o distinct -a $sample.scan.circRNA.psl.circRNA-exons.bed -b New_exon_and_Gencode_exon.bed > $sample.scan.circRNA.psl.genomic-exons.bed

# Prepare for comparison
cat $sample.scan.circRNA.psl.genomic-exons.bed | awk 'OFS="\t"{print $1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$4}' | sed 's/~/\t/g'| awk 'OFS="\t"{print $14,$2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | sortBed > temp.genomic-exons.bed

cat $sample.scan.circRNA.psl.annot.combine.txt | sortBed > $sample.scan.circRNA.psl.annot.combine.sort.txt

cat $sample.circRNA_candidates.annotated.txt | awk 'NR>1,OFS="\t"{print $2,$3,$4,$1,$6,$7,$10}' | sortBed > circRNA_candidates.annotated.txt.forExonAnalysis.bed
bedtools map -nonamecheck -f 0.95 -F 0.95 -c 4,7 -o distinct -a $sample.scan.circRNA.psl.annot.combine.sort.txt -b circRNA_candidates.annotated.txt.forExonAnalysis.bed | awk 'OFS="\t"{print $4,$2,$3,$1,0,"+",$10}' | sed 's/,circ/\t/g' | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7}' | sortBed > reads.annot.temp.bed

bedtools map -nonamecheck -c 7 -o distinct -a temp.genomic-exons.bed -b reads.annot.temp.bed | awk 'OFS="\t"{print $4,$2,$3,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' | sortBed > $sample.scan.circRNA.psl.genomic-exons.annot.bed

rm reads.annot.temp.bed temp.genomic-exons.bed $sample.scan.circRNA.psl.split.merge.flank2.posExons.bed $sample.scan.circRNA.psl.split.merge.flank2.negExons.bed
rm $sample.scan.circRNA.psl.annot.combine.sort.txt




echo "Starting the list with number of used exons for each circRNA"
date
tmp_dir=$(mktemp -d -t job-XXXXXXXXXX)
cat $sample.scan.circRNA.psl.genomic-exons.annot.bed | awk 'OFS="\t"{print $13,$14,$15}' | sort -k 1,1 -T $tmp_dir > temp1
cat temp1 | sort -k 2,2 -T $tmp_dir > temp2
cat temp2 | sort -k 3,3 -T $tmp_dir | uniq -c > $sample.scan.circRNA.psl.genomic-exons.annot.uniq.bed
rm -rf $tmp_dir




## New major improvement in version 6.1:
input=$sample.scan.circRNA.psl.genomic-exons.annot.uniq.bed

touch circRNA_exon_usage.txt
rm circRNA_exon_usage.txt

mkdir exon_usage_data

cat $input | awk '{print $4}' | sort | uniq | awk '{if ($1 != ".") print $1}'> circRNA-list
echo
echo "Getting circRNA exon usage"
while IFS='' read -r circRNA || [[ -n "$circRNA" ]]; do
        ## Getting the exons the specific circRNA checked now. Both novel and annotated exons are included, but  novel only if they have 50 or more reads.

#       ### If running on a system with SLURM queueing system (This will take a long time):
#	## Run 50 read version for panel
#	sbatch --account mtv_AU_project core_circ_grep_50reads.sh $circRNA $input $sample
#	## Run 10 read version for non-panel
#	sbatch --account mtv_AU_project core_circ_grep_10reads.sh $circRNA $input $sample
#
#	### If no SLURM queueing system - Or this script is already a sbatch script:
#        #Note: This removes novel exons with fewer than 10 reads.
#        grep $circRNA $input | awk '{print $3}' | sed 's/,/\n/g' | sort | uniq | grep -v "^[123456789]read_novelExon" > temp_exon_list
#        echo
#        echo "circRNA: "$circRNA
#        echo "exons in circRNA"
#        cat temp_exon_list
#        echo

	#Note: This removes novel exons with fewer than 50 reads. Only for panel datasets. Otherwise use cut off 10
        grep $circRNA $input | awk '{print $3}' | sed 's/,/\n/g' | sort | uniq | grep -v "^[123456789]read_novelExon" | grep -v "^[123][1234567890]read_novelExon" | grep -v "^4[123456789]read_novelExon" > temp_exon_list
        #echo
        echo "circRNA: "$circRNA > exon_usage.log
        #echo "exons in circRNA"
        #cat temp_exon_list
        #echo


        while IFS='' read -r exon || [[ -n "$exon" ]]; do
                exon_hit=$(grep $circRNA $sample.scan.circRNA.psl.genomic-exons.annot.uniq.bed | awk '{print $1,$2}' | grep $exon | awk '{split($0,a,","); sum += a[1]} END {print sum}')
                circRNA_coverage=$(grep $circRNA $sample.scan.circRNA.psl.genomic-exons.annot.uniq.bed | awk '{print $1,$3}' | grep $exon | awk '{split($0,a,","); sum += a[1]} END {print sum}')
                printf "$circRNA\t$exon_hit\t$circRNA_coverage\t$exon" | awk 'OFS="\t"{if($3 != 0) {print $1,$2,$3,$2/$3,$4;}}' >> exon_usage_data/$circRNA.circRNA_exon_usage.txt 2>> exon_usage.log

        done < temp_exon_list



done < circRNA-list

echo


#cat circRNA_exon_usage.txt | awk '$2>9' | awk '$4<0.9' | awk '$4>0.1' > circRNA_alternative_exon_usage.txt
#wc -l circRNA_exon_usage.txt circRNA_alternative_exon_usage.txt
#
#mkdir DELETE
#mv core_circ_grep.sh-*.out DELETE
#mv *.circRNA_exon_usage.txt DELETE

### after novel exons and alternative usage
cat exon_usage_data/*.circRNA_exon_usage.txt | sort -k 5,5 | sort -k 1,1 | uniq > $sample.circRNA_exon_usage.txt
# Filter to only keep rows with 5 columns. These are the real hits:
cat $sample.circRNA_exon_usage.txt | awk 'NF==5{print}{}' > $sample.circRNA_exon_usage_filter.txt

# Get the alternatively used exons
cat $sample.circRNA_exon_usage.txt | awk '$2>9' | awk '$4<0.9' | awk '$4>0.1' > $sample.circRNA_alternative_exon_usage.txt
#wc -l $sample.circRNA_exon_usage.txt $sample.circRNA_alternative_exon_usage.txt

### Making _circ_circRNA_exon_usage_length_of_exons
printf "Internal circRNA IDs\tExon used\tExon covered by read\tUsage level\texon\tstart\tend\tlength\n" > $sample.circ_circRNA_exon_usage_length_of_exons.txt
# only keep "Exon used" of 10 or more reads on exon. Then remove exons with name ".". Then sort by "Internal circRNA IDs"
cat $sample.circRNA_exon_usage_filter.txt | awk '$2>9' | awk '{ if ( $5 != "." ) { print $0; } }' | sort -k 1,1 > $sample.circ_circRNA_exon_usage_length_of_exons.temp.txt
cat $sample.circ_circRNA_exon_usage_length_of_exons.temp.txt | awk '{print $5}' | sed 's/_chr/\tchr/g' | awk '{print $2}' > coordinate.temp
cat coordinate.temp | sed 's/:/\t/g' | sed 's/-/\t/g' | awk 'OFS="\t"{print $2, $3, $3-$2}' > start_end_size
paste $sample.circ_circRNA_exon_usage_length_of_exons.temp.txt start_end_size >> $sample.circ_circRNA_exon_usage_length_of_exons.txt
rm coordinate.temp start_end_size $sample.circ_circRNA_exon_usage_length_of_exons.temp.txt


echo "Done with novel exons and alternative usage"
#echo "... Doing extra stuff v2"
echo
date







### Extra stuff v2

## Remove empty columns in circRNA candicates file:
mv $sample.circRNA_candidates.annotated.txt OLD.$sample.circRNA_candidates.annotated.txt
grep -v [[:space:]]0[[:space:]][+-][[:space:]]\.[[:space:]]\.[[:space:]]\.[[:space:]]\.[[:space:]]\. OLD.$sample.circRNA_candidates.annotated.txt > $sample.circRNA_candidates.annotated.txt


#wc -l OLD.$sample.circRNA_candidates.annotated.txt $sample.circRNA_candidates.annotated.txt


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


#echo "Cryptic novel exon comparison to annotated exons?"
#
#echo "novel exon sequences"


echo "Map read coverage on intron sequences"
echo "Get introns in circRNAs"
echo "This part was changed in v2"
## Make intron file:
# Getting unique introns that are located in circRNA expression regions
cat $intron_ucsc | awk 'OFS="\t"{print $1,$2,$3,"intron",0,$6}' | sortBed | uniq > introns.uniq.bed
cat $sample.circRNA_candidates.annotated.txt | grep -v internal_circRNA_name | awk 'OFS="\t"{print $2,$3,$4,$1,$6,$7}' | sortBed | bedtools map -c 4 -o distinct -a introns.uniq.bed -b - -nonamecheck > $sample.introns.uniq.circ.bed

# Remove introns that overlap an exon from Gencode
bedtools subtract -a introns.uniq.bed -b exon_annotation.reformat.bed -nonamecheck | sortBed | uniq | awk 'OFS="\t"{print $1,$2,$3,$4"_"$1":"$2"-"$3,1,$6}' | sortBed > introns.uniq.exon_remove.bed

echo "Getting read coverage in introns"
## This was changed in v2 so the intron retention is now only searched in bsj spanning reads, not all reads as before
bedtools coverage -bed -split -a introns.uniq.exon_remove.bed -b $sample.scan.circRNA.sort.bam > $sample.intron.coverage
echo "extracting the introns with >50% of nt covered by reads"
cat $sample.intron.coverage | awk '$10>0.5' > $sample.intron.coverage.50pct
# CircRNAs with more than 1% intronic coverage, if harboring a >50% coverage intron this is shown in last column
printf "internal_circRNA_name\tchr\tstart\tend\tBSJ_reads\tmean_read_coverage\tmean_intron_coverage\tintron_coverage\tcoverage_high_intron\n" > $sample.circRNA_intron_coverage.1pct.50pct_intron.txt
cat $sample.circRNA_intron_coverage.1pct.txt | awk 'OFS="\t"{print $2,$3,$4,$1,$5,$6,$7,$8}' | sortBed | bedtools map -F 1.0 -c 10 -o max -a - -b $sample.intron.coverage.50pct | awk 'OFS="\t"{print $4,$1,$2,$3,$5,$6,$7,$8,$9}'  >> $sample.circRNA_intron_coverage.1pct.50pct_intron.txt


bedtools coverage -a introns.uniq.exon_remove.bed -b $sample.psl.bed > $sample.introns.uniq.exon_remove.coverage.bed
cat $sample.circRNA_candidates.annotated.txt | grep -v internal_circRNA_name | awk 'OFS="\t"{print $2,$3,$4,$1,$6,$7}' | sortBed | bedtools map -c 4 -o distinct -a $sample.introns.uniq.exon_remove.coverage.bed -b - | awk 'OFS="\t"{print $0}' > $sample.introns.uniq.exon_remove.coverage.circ.bed
cat $sample.introns.uniq.exon_remove.coverage.circ.bed | grep circ_ >  $sample.introns.uniq.exon_remove.coverage.onlyCirc.bed
# Mapping novel exons on introns
# This step fails on bedtools==2.30.0!
mapBed -s -F 1.0 -c 4 -o distinct_only -a $sample.introns.uniq.exon_remove.coverage.onlyCirc.bed -b $sample.novel.exons.2reads.bed > $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.bed

# List of all unique introns in circRNA regions:
cat $sample.circRNA_candidates.annotated.txt | grep -v internal_circRNA_name | awk 'OFS="\t"{print $2,$3,$4,$1,$6,$7}' | intersectBed -s -u -a introns.uniq.bed -b - > $sample.all_circRNA_introns.bed

# This is a workaround for the bedtools >=2.30.0 error "***** ERROR: illegal number "1.0000000". Exiting..."
# By adding a 13th column we make sure bedtools does not assumes this is a valid BED12 file
# Bedtool's internal check will otherwise fail due to the float number at column 10 which is supposed to be an integer
# see https://github.com/arq5x/bedtools2/issues/981 for details

# old bedtools command was: bedtools map -F 1.0 -c 10 -o max -a $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.bed -b $sample.intron.coverage.50pct > $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.intronCov.bed

awk 'OFS="\t" {$13="."; print $0}' $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.bed  > $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.fixed.bed

bedtools map -F 1.0 -c 10 -o max -a $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.fixed.bed -b $sample.intron.coverage.50pct > $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.intronCov.bed

echo
echo "Number of introns"
wc -l $intron_ucsc introns.uniq.bed introns.uniq.exon_remove.bed $sample.introns.uniq.exon_remove.coverage.bed $sample.introns.uniq.exon_remove.coverage.circ.bed $sample.introns.uniq.exon_remove.coverage.onlyCirc.bed $sample.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.bed $sample.all_circRNA_introns.bed

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




if [[ $6 == "no" ]]
then
        echo "Cleaning up temporary output files"
        rm -rf exon_usage_data
        rm OLD.*

        #Remove everything that doesn't match this
        ls -I '*flank2.allExon*bed' -I '*circRNA_candidates.annotated.txt' \
          -I '*novel.exons.2reads.filter.bed' -I '*novel.exons.2reads.phases.tab' \
          -I '*novel.cryptic.spliced.exons.txt' -I '*circ_circRNA_exon_usage_length_of_exons.txt' \
          -I '*introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.intronCov.bed' \
          -I '*Potential_multi-round_circRNA.fa' -I '*Potential_multi-round_circRNA.psl.annot*' \
          -I '*bam' -I '*bai' -I '*bw' | xargs rm

else
        echo "Keeping all of the temporary output files"
fi


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


