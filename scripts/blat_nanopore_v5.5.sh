#!/bin/bash
#SBATCH -c 8
#SBATCH -p normal
#SBATCH --mem=128g
#SBATCH --time=72:00:00

# This is version 5.5 of the Nanopore circRNA dataction pipeline
# This version uses blat to map reads.  NanoFilt is used to trim low quality reads (q 7). Parallel BLAT is used to map nanopore reads to genome.



data_folder=$1
sample=$2
species=$3
reference_path=$4
scriptFolder=$5

mkdir $6
cd $6

#OriginalScriptFolder=~/long_read_circRNA_v2/scripts
#mkdir scripts

#cp -r $OriginalScriptFolder/* $scriptFolder


if [ "$species" = "human" ] ; then
fa=$reference_path/human/hg19.fa
mRNA=$reference_path/human/Human_refFlat_hg19_Oct2018.unique.merge.bed
exon=$reference_path/human/Human_refFlat_exon_hg19_Oct2018.merge.bed
single_exon=$reference_path/human/Human_refFlat_exon_hg19_Oct2018.sort.bed
est=$reference_path/human/UCSC-EST-exons_hg19_09-2018.bed
circBase=$reference_path/human/hsa_circRNA_complete.hg19.unique.sort.length.bed
circAtlas=$reference_path/human/circAtlas2.0_June2019_human_hg19_circRNA.0-based.bed
CIRCpedia=$reference_path/human/CIRCpedia_v2_June2019_human_hg19_All_circRNA.unique.bed
genomeSize=$reference_path/human/hg19.chrom.sizes
circRNA_prefix=hsa_circ_
fi

if [ "$species" = "mouse" ] ; then
fa=$reference_path/mouse/Mus_musculus.GRCm38.87.chr-fix.fa
mRNA=$reference_path/mouse/Mouse_refFlat_mm10_Oct2018.unique.merge.bed
exon=$reference_path/mouse/Mouse_refFlat_exon_mm10_Oct2018.merge.bed
single_exon=$reference_path/mouse/Mouse_refFlat_exon_mm10_Oct2018.sort.bed
est=$reference_path/mouse/UCSC-EST-exons_mm10_09-2018.bed
circBase=$reference_path/mouse/mmu_circRNA_complete.mm10lift.sort.length.bed
circAtlas=$reference_path/mouse/circAtlas2.0_Aug2019_mouse_mm10_circRNA.0-based.bed
CIRCpedia=$reference_path/mouse/CIRCpedia_v2_June2019_mouse_mm10_All_circRNA.unique.bed
genomeSize=$reference_path/mouse/mm10.chrom.sizes
circRNA_prefix=mmu_circ_
fi


# Temp folder
temp_sort=$(mktemp -d /tmp/foo.XXXXXXXXX)

date

echo "Sample: "$sample
echo
echo "Number of raw reads before NanoFilt -q 7 -l 250"
zcat $data_folder/$sample.fq.gz | wc -l | awk '{print $1/4}'
#echo
#echo "Dates when reads were sequenced:"
#zcat $data_folder/$sample.fq.gz | grep "start_time=" | sed 's/start_time=/;/g' | cut -d \; -f 2 | sed 's/T/ /g'  | awk '{print $1}'  | sort | uniq -c | sort -nrk 1,1

echo
date
echo "NanoFilt to remove reads under quality 7 and conversion to fasta"
zcat $data_folder/$sample.fq.gz | NanoFilt -q 7 -l 250 | sed -n '1~4s/^@/>/p;2~4p' > $sample.fa
echo
date
echo "Number of filtered reads after NanoFilt -q 7 -l 250"
cat $sample.fa | wc -l | awk '{print $1/2}'
echo
date
echo "Mapping with pblat - parallelized blat with multi-threads support (http://icebert.github.io/pblat/)"
echo "lower case sequences in the genome file are masked out"
echo "Showing a dot for every 50k sequences processed"
pblat  -threads=8 -trimT -dots=50000 -mask=lower $fa $sample.fa $sample.psl
echo "Blat done"
date

# Remove non-standard chromosomes
cat $sample.psl | grep -v "_random" | grep -v "_hap" | grep -v "chrUn_" > $sample.temp.psl
rm $sample.psl
mv $sample.temp.psl $sample.psl
cat $sample.psl | awk '{print $10}' | sort | uniq -c | sort -nrk 1,1 > mappings_per_read.txt
cat mappings_per_read.txt | awk '{print $1}' | uniq -c | sort -nrk 1,1 > $sample.histogram_number_of_genomic_hits_per_read.txt
cat $sample.psl | perl $scriptFolder/psl2bed12.pl | sortBed > $sample.psl.bed


echo
date
echo "Getting group numbers"
cat $sample.psl | perl $scriptFolder/blat_output_processing_v2.pl > $sample.scan.psl

# Read fragment numbers
echo
date
echo "The different groups, numbers of read fragments:"
cat $sample.scan.psl | awk '{print $NF}' | sort | uniq -c | sort -nrk 1,1 | head -6
cat $sample.scan.psl | awk '{print $NF}' | sort | uniq -c | sort -nrk 1,1 | head -6 > $sample.scan.groupNumbers.fragments.txt
# Read numbers
echo
echo "The different groups, numbers of unique reads:"
cat $sample.scan.psl | awk '{print $10,$NF}' | sort | uniq | awk '{print $NF}' | sort | uniq -c | sort -nrk 1,1 | head -6
cat $sample.scan.psl | awk '{print $10,$NF}' | sort | uniq | awk '{print $NF}' | sort | uniq -c | sort -nrk 1,1 | head -6 > $sample.scan.groupNumbers.reads.txt


# Get circRNA reads
head -5 $sample.scan.psl > $sample.scan.circRNA.psl
grep circRNA $sample.scan.psl >> $sample.scan.circRNA.psl
## converting psl to bed12
cat $sample.scan.circRNA.psl | perl $scriptFolder/psl2bed12.pl | sortBed > $sample.scan.circRNA.psl.bed



echo
echo "Making bam and bigWig (.bw) files for use in genome browsers"
# All Blat mapped reads
bedtools bedtobam -i $sample.psl.bed -bed12 -g $genomeSize > $sample.bam
samtools sort $sample.bam > $sample.sort.bam
samtools index $sample.sort.bam
rm $sample.bam
bamCoverage --binSize 1 --numberOfProcessors 8 -b $sample.sort.bam -o $sample.bw

# Blat mapped BSJ spanning reads
bedtools bedtobam -i $sample.scan.circRNA.psl.bed -bed12 -g $genomeSize > $sample.circRNA.bam
samtools sort $sample.circRNA.bam > $sample.circRNA.sort.bam
samtools index $sample.circRNA.sort.bam
rm $sample.circRNA.bam
bamCoverage --binSize 1 --numberOfProcessors 8 -b $sample.circRNA.sort.bam -o $sample.circRNA.bw



echo
date
echo "outputting Potential_multi-round_circRNA"
# Get Potential_multi-round_circRNA
head -5 $sample.scan.psl > $sample.scan.Potential_multi-round_circRNA.psl
grep Potential_multi-round_circRNA $sample.scan.psl >> $sample.scan.Potential_multi-round_circRNA.psl
## converting psl to bed12
cat $sample.scan.Potential_multi-round_circRNA.psl | perl $scriptFolder/psl2bed12.pl | sortBed > $sample.scan.Potential_multi-round_circRNA.psl.bed
bedtools bedtobam -i $sample.scan.Potential_multi-round_circRNA.psl.bed -bed12 -g $genomeSize > $sample.scan.Potential_multi-round_circRNA.bam
samtools sort $sample.scan.Potential_multi-round_circRNA.bam > $sample.scan.Potential_multi-round_circRNA.sort.bam
samtools index $sample.scan.Potential_multi-round_circRNA.sort.bam
# Making a special file to show how many rounds each read takes
cat $sample.scan.Potential_multi-round_circRNA.psl.bed | sed 's/~/\t/g' | awk 'OFS="\t"{print $4,$2,$3,$1,$6,$7}' | sortBed | awk 'OFS="\t"{print $1"~"$4,$2,$3,$4,$5,$6}' | bedtools merge > $sample.scan.Potential_multi-round_circRNA.psl.merge.bed
cat $sample.scan.Potential_multi-round_circRNA.psl.bed | sed 's/~/\t/g' | awk 'OFS="\t"{print $4,$2,$3,$1,$6,$7}' | sortBed | awk 'OFS="\t"{print $1"~"$4,$2,$3,$4,$5,$6}' | bedtools coverage -counts -a $sample.scan.Potential_multi-round_circRNA.psl.merge.bed -b - > temp.$sample.multi-round.count.txt
printf "#Chr\tStart\tEnd\tRead_name\tNumber_of_rounds\tOverlapping_gene\n" > $sample.scan.Potential_multi-round_circRNA.psl.annot.bed
cat temp.$sample.multi-round.count.txt | sed 's/~/\t/g' | awk 'OFS="\t"{print $2,$3,$4,$1,$5}' | sortBed | bedtools map -c 4 -o distinct -a - -b $mRNA >> $sample.scan.Potential_multi-round_circRNA.psl.annot.bed
cat $sample.scan.Potential_multi-round_circRNA.psl.annot.bed | awk '{print $NF}' | sort | uniq -c | sort -nrk 1,1 > $sample.scan.Potential_multi-round_circRNA.psl.annot.count.txt
cat $sample.scan.Potential_multi-round_circRNA.psl.annot.bed | awk '{print $4}' | grep -v Read_name > $sample.temp.read_names
#grep --no-group-separator -A1 -f $sample.temp.read_names $sample.fa > $sample.Potential_multi-round_circRNA.fa
samtools faidx $sample.fa
echo "Using samtools to make Potential_multi-round_circRNA fasta files for sample: "$sample.fa
xargs samtools faidx $sample.fa < $sample.temp.read_names > $sample.Potential_multi-round_circRNA.fa
rm temp.$sample.multi-round.count.txt $sample.temp.read_names

echo
date
echo "Annotating circRNAs and converting to bed6 format"

### Adding split
echo "  First splitting bed file in 8 to optimize run time"
bedtools split -n 8 -p $sample.scan.circRNA.psl -a simple -i $sample.scan.circRNA.psl.bed

## Very time consuming step for large datasets... Converting to bed6 (turning each block from bed12 into a single bed entry). Header is removed using awk
cat $sample.scan.circRNA.psl.00001.bed | bed12ToBed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00001.bed &
cat $sample.scan.circRNA.psl.00002.bed | bed12ToBed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00002.bed &
cat $sample.scan.circRNA.psl.00003.bed | bed12ToBed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00003.bed &
cat $sample.scan.circRNA.psl.00004.bed | bed12ToBed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00004.bed &
cat $sample.scan.circRNA.psl.00005.bed | bed12ToBed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00005.bed &
cat $sample.scan.circRNA.psl.00006.bed | bed12ToBed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00006.bed &
cat $sample.scan.circRNA.psl.00007.bed | bed12ToBed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00007.bed &
cat $sample.scan.circRNA.psl.00008.bed | bed12ToBed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00008.bed &

wait

cat $sample.scan.circRNA.psl.annot.00001.bed $sample.scan.circRNA.psl.annot.00002.bed $sample.scan.circRNA.psl.annot.00003.bed $sample.scan.circRNA.psl.annot.00004.bed $sample.scan.circRNA.psl.annot.00005.bed $sample.scan.circRNA.psl.annot.00006.bed $sample.scan.circRNA.psl.annot.00007.bed $sample.scan.circRNA.psl.annot.00008.bed > $sample.scan.circRNA.psl.annot.bed
rm $sample.scan.circRNA.psl.annot.0*.bed

### Done with split



printf "#chr\tstart\tend\tread_name\tread_length\tgene_coverage\texon_coverage\tEST_coverage\tintron_coverage\n" > $sample.scan.circRNA.psl.annot.txt
cat $sample.scan.circRNA.psl.annot.bed | sort -k 4,4 -T $temp_sort | awk 'OFS="\t"{print $1,$2,$3,$4,$3-$2,($3-$2)*$8,($3-$2)*$10,($3-$2)*$12,($3-$2)*$8-($3-$2)*$10}' | perl $scriptFolder/combine_annot_segments.pl >> $sample.scan.circRNA.psl.annot.txt
cat $sample.scan.circRNA.psl.annot.txt | perl $scriptFolder/make_circRNAs_from_annot.txt.pl > $sample.scan.circRNA.psl.annot.combine.txt

echo
date
echo "Refining circRNA edges based annotated exon boundaries and annotated circRNAs"
# Making a unique list of circRNAs
cat $sample.scan.circRNA.psl.annot.combine.txt | awk 'NR>1,OFS="\t"{print $1,$2,$2+1,$4"~"$5"~"$6"~"$7"~"$8"~"$9}' | sortBed | uniq > temp_start
cat $sample.scan.circRNA.psl.annot.combine.txt | awk 'NR>1,OFS="\t"{print $1,$3-1,$3,$4"~"$5"~"$6"~"$7"~"$8"~"$9}' | sortBed | uniq > temp_end
#sortBed -i $single_exon > exon_ref
# Prints the start and end position of closest exon. In special cases where the circRNA is produced far inside an annoteted exon, such as occurs for Malat1, are filtered away.
bedtools closest -t first -d -header -a temp_start -b $single_exon -nonamecheck | awk 'OFS="\t"{if ($2 - $6 < 31) print $1,$6,$7,$4,$10,$11}' | awk 'OFS="\t"{if($2 > -1) print $0 }' > temp_start.exon
bedtools closest -t first -d -header -a temp_end -b $single_exon -nonamecheck | awk 'OFS="\t"{if ($7 - $3 < 31) print $1,$6,$7,$4,$10,$11}' | awk 'OFS="\t"{if($2 > -1) print $0 }' > temp_end.exon
cat temp_start.exon temp_end.exon | sort -k 4,4 > temp_edge_exon
bedtools groupby -g 4 -c 1,2,3,5,6,4 -o distinct,min,max,distinct,max,count -i temp_edge_exon | awk 'OFS="\t"{print $2,$3,$4,$1,$6,$5,$7}' > temp_exon-ends0
# Allowing only edges that are formed from 2 read segments:
cat temp_exon-ends0 | awk 'OFS="\t"{if ($7 == 2) print $1,$2,$3,$4,$5,$6}' > temp_exon-ends
## Exon match: Max 30 bp distance, correct strand
grep -v "+,-" temp_exon-ends | awk 'OFS="\t"{if($5 < 31) print $0 }' > $sample.scan.circRNA.psl.annot.combine.correct.bed
cat $sample.scan.circRNA.psl.annot.combine.correct.bed | awk 'OFS="\t"{print $1,$2,$3,"exon_match",0,$6}' | sort -nk 3,3 | sort -nk 2,2 | sort -k 1,1 | uniq  | sortBed > base_list_exon-match.bed
cat $sample.scan.circRNA.psl.annot.combine.correct.bed | awk 'OFS="\t"{print $1,$2,$3,$5,$6,$4}' | sed 's/~/\t/g' | awk 'OFS="\t"{print $1,$2,$3,$6,$4,$5,$7,$8,$9,$10,$11}' | sortBed > $sample.scan.circRNA.psl.annot.combine.correct.full.bed
## No exon match: Over 30 bp distance or segments on different strands
grep "+,-" temp_exon-ends | awk 'OFS="\t"{ print $4 }' > temp_exon-ends_nohit
cat temp_exon-ends | awk 'OFS="\t"{if($5 > 30) print $4 }' >> temp_exon-ends_nohit
cat temp_exon-ends_nohit | sed 's/~/\t/g' | awk 'OFS="\t"{print $1}' | sort | uniq > temp_exon-ends_nohit_uniq

#rm temp_start temp_end exon_ref temp_start.exon temp_end.exon temp_edge_exon temp_exon-ends temp_exon-ends_nohit temp_exon-ends0

### For the base_list_exon-match.bed file
# finding host gene: Any overlap of the circRNA with an annotated refSeq gene on the SENSE strand
bedtools map -s -c 4 -o distinct -a base_list_exon-match.bed -b $mRNA -nonamecheck > base_list_exon-match.temp
#intersectBed -wao -f 0.5 -s -a base_list_exon-match.bed -b $mRNA | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{print $NF,$0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{print $0,$1}' | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' | sed s'/\t\t/\t/g' > base_list_exon-match.temp
# finding antisense genes: ANY overlap of circRNA with an annotated refSeq gene on the ANTISENSE strand
#intersectBed -wao -S -a base_list_exon-match.temp -b $mRNA | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{print $NF,$0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{print $0,$1}' | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' | sed 's/\t\t/\t/'g > base_list_exon-match.temp2
# Mapping circBase ID
bedtools map -f 1.0 -F 1.0 -c 4 -o distinct -a base_list_exon-match.temp -b $circBase -nonamecheck | awk '!($0 in a) {a[$0];print}' > base_list_exon-match.temp2.bed
# Mapping circAtlas
bedtools map -f 1.0 -F 1.0 -c 4 -o distinct -a base_list_exon-match.temp2.bed -b $circAtlas -nonamecheck | awk '!($0 in a) {a[$0];print}' > base_list_exon-match.temp3.bed
# Mapping CIRCpedia
bedtools map -f 1.0 -F 1.0 -c 4 -o distinct -a base_list_exon-match.temp3.bed -b $CIRCpedia -nonamecheck | awk '!($0 in a) {a[$0];print}' > base_list_exon-match.temp4.bed

# Mapping stuff on to the base list
bedtools map -f 1.0 -F 1.0 -c 4,7,8,9,10,11,5,5,5 -o count,mean,mean,mean,mean,mean,min,max,mean -a base_list_exon-match.temp4.bed -b $sample.scan.circRNA.psl.annot.combine.correct.full.bed -nonamecheck | awk 'OFS="\t"{print $1,$2,$3,$4,$11,$6,$7,$8,$9,$10,$12,$13,$14,$15,$16,$17,$18,$19}' | sort -nrk 5,5 > base_list_exon-match.annot.prefilter.bed
echo
#echo "Removing exon_match in chrM base_list_exon-match"
#grep chrM base_list_exon-match.annot.prefilter.bed | wc -l
#echo "Removing exon_match to Rn45s from base_list_exon-match"
#grep Rn45s base_list_exon-match.annot.prefilter.bed | wc -l
cat base_list_exon-match.annot.prefilter.bed | grep -v chrM | grep -v Rn45s > $sample.base_list_exon-match.annot.bed


rm base_list_exon-match.temp base_list_exon-match.temp3.bed base_list_exon-match.temp4.bed base_list_exon-match.annot.prefilter.bed


# For v 5.5 I increased the stringency. See below
### for the reads that do not match exons I check for similarity to circBase, circAtlas or CIRCpedia circRNA. If this is found, the annotated circRNA entry defines boundaries.
cat $sample.scan.circRNA.psl.annot.combine.txt | sortBed > $sample.scan.circRNA.psl.annot.combine.sort.txt
# finding host gene: Any overlap of the circRNA with an annotated refSeq gene on either
intersectBed -wao -a $sample.scan.circRNA.psl.annot.combine.sort.txt -b $mRNA -nonamecheck | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{print $NF,$0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{print $0,$1}' | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' | sed s'/\t\t/\t/g' > $sample.scan.circRNA.psl.annot.combine.sort.temp


## FIX in v 5.5! (was potentially misannotation long circRNAs like RMST)
# No longer Annotating with 95% similar circRNAs
# Now we first annotate with 99.9%, then 99% then 95% similar circRNAs. The highest similarity that generates a hit wins. If no hit for 95% then there is just a dot
bedtools map -f 0.999 -F 0.999 -c 4 -o last -a $sample.scan.circRNA.psl.annot.combine.sort.temp -b $circBase -nonamecheck | awk '!($0 in a) {a[$0];print}' > $sample.scan.circRNA.psl.annot.combine.sort.temp2a
bedtools map -f 0.99 -F 0.99 -c 4 -o last -a $sample.scan.circRNA.psl.annot.combine.sort.temp2a -b $circBase -nonamecheck | awk '!($0 in a) {a[$0];print}' > $sample.scan.circRNA.psl.annot.combine.sort.temp2b
bedtools map -f 0.95 -F 0.95 -c 4 -o last -a $sample.scan.circRNA.psl.annot.combine.sort.temp2b -b $circBase -nonamecheck | awk '!($0 in a) {a[$0];print}' > $sample.scan.circRNA.psl.annot.combine.sort.temp2c
cat $sample.scan.circRNA.psl.annot.combine.sort.temp2c | awk 'OFS="\t"{print $0, $(NF-2)"space"$(NF-1)"space"$NF}' | sed "s/\.space//g" | sed "s/space.*//g" | awk 'OFS="\t"{$(NF-3)=$(NF-2)=$(NF-1)="";print $0}' | sed 's/\t\+/\t/g;s/^\t//' > $sample.scan.circRNA.psl.annot.combine.sort.temp2
#rm $sample.scan.circRNA.psl.annot.combine.sort.temp2a $sample.scan.circRNA.psl.annot.combine.sort.temp2b $sample.scan.circRNA.psl.annot.combine.sort.temp2c
bedtools map -f 0.999 -F 0.999 -c 4 -o last -a $sample.scan.circRNA.psl.annot.combine.sort.temp2 -b $circAtlas -nonamecheck | awk '!($0 in a) {a[$0];print}' > $sample.scan.circRNA.psl.annot.combine.sort.temp3a
bedtools map -f 0.99 -F 0.99 -c 4 -o last -a $sample.scan.circRNA.psl.annot.combine.sort.temp3a -b $circAtlas -nonamecheck | awk '!($0 in a) {a[$0];print}' > $sample.scan.circRNA.psl.annot.combine.sort.temp3b
bedtools map -f 0.95 -F 0.95 -c 4 -o last -a $sample.scan.circRNA.psl.annot.combine.sort.temp3b -b $circAtlas -nonamecheck | awk '!($0 in a) {a[$0];print}' > $sample.scan.circRNA.psl.annot.combine.sort.temp3c
cat $sample.scan.circRNA.psl.annot.combine.sort.temp3c | awk 'OFS="\t"{print $0, $(NF-2)"space"$(NF-1)"space"$NF}' | sed "s/\.space//g" | sed "s/space.*//g" | awk 'OFS="\t"{$(NF-3)=$(NF-2)=$(NF-1)="";print $0}' | sed 's/\t\+/\t/g;s/^\t//' > $sample.scan.circRNA.psl.annot.combine.sort.temp3
#rm $sample.scan.circRNA.psl.annot.combine.sort.temp3a $sample.scan.circRNA.psl.annot.combine.sort.temp3b $sample.scan.circRNA.psl.annot.combine.sort.temp3c
bedtools map -f 0.999 -F 0.999 -c 4 -o last -a $sample.scan.circRNA.psl.annot.combine.sort.temp3 -b $CIRCpedia -nonamecheck | awk '!($0 in a) {a[$0];print}' > $sample.scan.circRNA.psl.annot.combine.circID.beda
bedtools map -f 0.99 -F 0.99 -c 4 -o last -a $sample.scan.circRNA.psl.annot.combine.circID.beda -b $CIRCpedia -nonamecheck | awk '!($0 in a) {a[$0];print}' > $sample.scan.circRNA.psl.annot.combine.circID.bedb
bedtools map -f 0.95 -F 0.95 -c 4 -o last -a $sample.scan.circRNA.psl.annot.combine.circID.bedb -b $CIRCpedia -nonamecheck | awk '!($0 in a) {a[$0];print}' > $sample.scan.circRNA.psl.annot.combine.circID.bedc
cat $sample.scan.circRNA.psl.annot.combine.circID.bedc | awk 'OFS="\t"{print $0, $(NF-2)"space"$(NF-1)"space"$NF}' | sed "s/\.space//g" | sed "s/space.*//g" | awk 'OFS="\t"{$(NF-3)=$(NF-2)=$(NF-1)="";print $0}' | sed 's/\t\+/\t/g;s/^\t//' > $sample.scan.circRNA.psl.annot.combine.circID.bed
#rm $sample.scan.circRNA.psl.annot.combine.circID.beda $sample.scan.circRNA.psl.annot.combine.circID.bedb $sample.scan.circRNA.psl.annot.combine.circID.bedc

# Getting only the no_match reads:
grep -Fwf temp_exon-ends_nohit_uniq $sample.scan.circRNA.psl.annot.combine.circID.bed | sortBed > no_exon_match_reads.bed

#circBase
cat no_exon_match_reads.bed | awk 'OFS="\t"{print $11,$12,$13}' | uniq | grep -v "^\." | awk '{print $1}' | sed 's/,/\t/g' | awk 'OFS="\t"{print $1}' | sort | uniq | grep -v "\." > circBase.no_exon_match_reads_circIDs.txt
#circAtlas
cat no_exon_match_reads.bed | awk 'OFS="\t"{print $11,$12,$13}' | uniq | grep "^\."| grep -v "^\.[[:space:]]\." | awk '{print $2}' | sed 's/,/\t/g' | awk 'OFS="\t"{print $1}' | sort | uniq | grep -v "\." > circAtlas.no_exon_match_reads_circIDs.txt
#CIRCpedia
cat no_exon_match_reads.bed | awk 'OFS="\t"{print $11,$12,$13}' | uniq | grep "^\."| grep "^\.[[:space:]]\." | grep -v "^\.[[:space:]]\.[[:space:]]\." | awk '{print $3}' | sed 's/,/\t/g' | awk 'OFS="\t"{print $1}' | sort | uniq | grep -v "\." > CIRCpedia.no_exon_match_reads_circIDs.txt
# Getting the genome coordinates from circBase from the no-exon aligning reads
grep -Fwf circBase.no_exon_match_reads_circIDs.txt $circBase > $sample.circBase_no_exon_match_reads.bed
# Getting the genome coordinates from circAtlas from the no-exon aligning reads
grep -Fwf circAtlas.no_exon_match_reads_circIDs.txt $circAtlas >> $sample.circBase_no_exon_match_reads.bed
# Getting the genome coordinates from CIRCpedia from the no-exon aligning reads
grep -Fwf CIRCpedia.no_exon_match_reads_circIDs.txt $circBase >> $sample.circBase_no_exon_match_reads.bed
cat $sample.circBase_no_exon_match_reads.bed | sortBed > temp
mv temp $sample.circBase_no_exon_match_reads.bed

# Mapping stuff on to the $sample.circBase_no_exon_match_reads.bed list - I increases the stringency from 95 to 99% here for v 5.1:
bedtools map -f 0.99 -F 0.99 -c 4,10,11,12,13,5,6,7,8,9 -o count,distinct,distinct,distinct,distinct,mean,mean,mean,mean,mean -a $sample.circBase_no_exon_match_reads.bed -b no_exon_match_reads.bed -nonamecheck | awk 'OFS="\t"{print $1,$2,$3,"no_exon_match",$7,$6,$8,$9,$10,$11,$12,$13,$14,$15,$16,"NA","NA","NA"}' | sort -nrk 5,5 > base_list_no-exon.cirBaseID.annot.prefilter.bed
cat base_list_no-exon.cirBaseID.annot.prefilter.bed | grep -v chrM | grep -v Rn45s | grep -v "no_exon_match[[:space:]]0" > $sample.base_list_no-exon.cirBaseID.annot.bed

#rm base_list_no-exon.cirBaseID.annot.prefilter.bed

echo
date
echo "Generating internal_circRNA_name and outputting candidate circRNA list"
# Combine the positive hits
cat $sample.base_list_exon-match.annot.bed $sample.base_list_no-exon.cirBaseID.annot.bed | sort -nrk 5,5 > temp.circ.hits
printf "chr\tstart\tend\tdescription\tBSJ_reads\tstrand\tgene\tcircBase_ID\tcircAtlas_ID\tCIRCpedia_ID\tmean_read_coverage\tmean_gene_coverage\tmean_exon_coverage\tmean_EST_coverage\tmean_intron_coverage\tmin_exon_adjust\tmax_exon_adjust\tmean_exon_adjust\n" > $sample.circRNA_candidates.annotated.bed
cat temp.circ.hits >> $sample.circRNA_candidates.annotated.bed

        count=0
        while read -r line
        do
                       if [ $count == 0 ]; then
                               echo "internal_circRNA_name" > circRNA_name.temp
                       else
                               name=$(echo $line | awk '{print $7}')
                               #antisense_name=$(echo $line | awk '{print $8}')
                               if [ $name = "." ] || [ $name = "NA" ]; then
                               #        if [ $antisense_name = "." ] || [ $antisense_name = "NA" ]; then
                                               circ_host="intergenic"
                               #        else
                               #                circ_host=$antisense_name"-AS"
                               #        fi
                               else
                                       circ_host=$name
                               fi
                               echo `printf $circRNA_prefix%04d $count`"_$circ_host" >> circRNA_name.temp
                       fi
               count=$(expr $count + 1)
        done < $sample.circRNA_candidates.annotated.bed
	paste circRNA_name.temp $sample.circRNA_candidates.annotated.bed > $sample.circRNA_candidates.annotated.txt


#### for the reads that do not match exons and also does not have 99% similarity to known circRNAs
cat no_exon_match_reads.bed | uniq | grep "\.[[:space:]]\.[[:space:]]\." | sortBed > no_exon_no_circRNA.bed


## Delete temp files
rm temp.circ.hits $sample.scan.circRNA.psl.annot.combine.sort.temp temp_exon-ends_nohit_uniq $sample.scan.circRNA.psl.annot.combine.correct.full.bed circRNA_name.temp $sample.scan.circRNA.psl.annot.combine.correct.bed  $sample.scan.circRNA.psl.annot.combine.sort.temp2 $sample.scan.circRNA.psl.annot.combine.sort.temp3
rm -r $temp_sort
rm temp* $sample.scan.Potential_multi-round_circRNA.bam $sample.scan.Potential_multi-round_circRNA.psl.bed $sample.scan.Potential_multi-round_circRNA.psl
#rm $sample.scan.circRNA.psl.annot.combine.circID.bed $sample.scan.circRNA.psl.annot.combine.sort.txt $sample.base_list_exon-match.annot.bed
rm $sample.scan.Potential_multi-round_circRNA.psl.merge.bed $sample.scan.Potential_multi-round_circRNA.sort.bam.bai
rm $sample.scan.Potential_multi-round_circRNA.sort.bam
#rm $sample.scan.circRNA.bam $sample.scan.circRNA.psl.annot.bed $sample.scan.circRNA.psl.annot.txt $sample.scan.circRNA.psl $sample.sort.bam $sample.scan.circRNA.psl.bed
#rm $sample.sort.bam.bai $sample.psl.bed $sample.fa.fai $sample.psl $sample.fa no_exon_no_circRNA.bed
rm base_list_exon-match.bed no_exon_match_reads.bed circBase.no_exon_match_reads_circIDs.txt circAtlas.no_exon_match_reads_circIDs.txt CIRCpedia.no_exon_match_reads_circIDs.txt
#rm base_list_no-exon.cirBaseID.annot.prefilter.bed
rm mappings_per_read.txt $sample.scan.psl $sample.circBase_no_exon_match_reads.bed #$sample.base_list_no-exon.cirBaseID.annot.bed

echo
date
echo
echo
echo Done

