#!/bin/bash
#SBATCH -c 1
#SBATCH -p normal
#SBATCH --mem=32g
#SBATCH --time=12:00:00


species=human
sample=20191230_20200102_hBrain_mRNA

if [ "$species" = "human" ] ; then
fa=~/long_read_circRNA_v2/data/human/hg19.fa
mRNA=~/long_read_circRNA_v2/data/human/Human_refFlat_hg19_Oct2018.unique.merge.bed
exon=~/long_read_circRNA_v2/data/human/Human_refFlat_exon_hg19_Oct2018.merge.bed
single_exon=~/long_read_circRNA_v2/data/human/Human_refFlat_exon_hg19_Oct2018.sort.bed
est=~/long_read_circRNA_v2/data/human/UCSC-EST-exons_hg19_09-2018.bed
circBase=~/long_read_circRNA_v2/data/human/hsa_circRNA_complete.hg19.unique.sort.length.bed
circAtlas=~/long_read_circRNA_v2/data/human/circAtlas2.0_June2019_human_hg19_circRNA.0-based.bed
CIRCpedia=~/long_read_circRNA_v2/data/human/CIRCpedia_v2_June2019_human_hg19_All_circRNA.unique.bed
genomeSize=~/long_read_circRNA_v2/data/human/hg19.chrom.sizes
circRNA_prefix=hsa_circ_
fi

if [ "$species" = "mouse" ] ; then
fa=~/long_read_circRNA_v2/data/mouse/Mus_musculus.GRCm38.87.chr-fix.fa
mRNA=~/long_read_circRNA_v2/data/mouse/Mouse_refFlat_mm10_Oct2018.unique.merge.bed
exon=~/long_read_circRNA_v2/data/mouse/Mouse_refFlat_exon_mm10_Oct2018.merge.bed
single_exon=~/long_read_circRNA_v2/data/mouse/Mouse_refFlat_exon_mm10_Oct2018.sort.bed
est=~/long_read_circRNA_v2/data/mouse/UCSC-EST-exons_mm10_09-2018.bed
circBase=~/long_read_circRNA_v2/data/mouse/mmu_circRNA_complete.mm10lift.sort.length.bed
circAtlas=~/long_read_circRNA_v2/data/mouse/circAtlas2.0_Aug2019_mouse_mm10_circRNA.0-based.bed
CIRCpedia=~/long_read_circRNA_v2/data/mouse/CIRCpedia_v2_June2019_mouse_mm10_All_circRNA.unique.bed
genomeSize=~/long_read_circRNA_v2/data/mouse/mm10.chrom.sizes
circRNA_prefix=mmu_circ_
fi

scriptFolder=~/long_read_circRNA_v2/scripts


#date
#echo "make1col.fa"
#cat $sample.fa | awk '{print $1}' > $sample.1col.fa

#date
#echo "make psl.bed"
#cat $sample.scan.circRNA.psl | perl $scriptFolder/psl2bed12.pl > $sample.scan.circRNA.psl.bed
echo "make scan.circRNA.psl.query.bed"
cat $sample.scan.circRNA.psl | perl $scriptFolder/psl2bed12_query_v2.pl > $sample.scan.circRNA.psl.query.bed
date
echo "make scan.circRNA.psl.query.seq.tab"
bedtools getfasta -s -tab -split -nameOnly -fi $sample.1col.fa -bed $sample.scan.circRNA.psl.query.bed > $sample.scan.circRNA.psl.query.seq.tab
date
echo "make scan.circRNA.psl.bam"
bedtools bedtobam -i $sample.scan.circRNA.psl.bed -bed12 -g $genomeSize > $sample.scan.circRNA.psl.bam
date
echo "make scan.circRNA.psl.sam"
samtools view -H $sample.scan.circRNA.psl.bam > $sample.scan.circRNA.psl.sam.header
samtools view $sample.scan.circRNA.psl.bam > $sample.scan.circRNA.psl.sam
date
echo "make scan.circRNA.psl.seq.bam"
paste $sample.scan.circRNA.psl.sam $sample.scan.circRNA.psl.query.seq.tab | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13,$11}' > $sample.scan.circRNA.psl.query.temp.seq.sam
cat $sample.scan.circRNA.psl.sam.header $sample.scan.circRNA.psl.query.temp.seq.sam > $sample.scan.circRNA.psl.query.temp.seq.header.sam
samtools view -b $sample.scan.circRNA.psl.query.temp.seq.header.sam > $sample.scan.circRNA.psl.seq.bam
date
echo "sort and index"
samtools sort -O bam -T temp.$sample $sample.scan.circRNA.psl.seq.bam > $sample.scan.circRNA.psl.seq.sort.bam
samtools index $sample.scan.circRNA.psl.seq.sort.bam


#date
#echo "make psl.bed"
#cat $sample.psl | perl $scriptFolder/psl2bed12.pl > $sample.psl.bed
echo "make psl.query.bed"
cat $sample.psl | perl $scriptFolder/psl2bed12_query_v2.pl > $sample.psl.query.bed
date
echo "make psl.query.seq.tab"
bedtools getfasta -s -tab -split -nameOnly -fi $sample.1col.fa -bed $sample.psl.query.bed > $sample.psl.query.seq.tab
date
echo "make psl.bam"
bedtools bedtobam -i $sample.psl.bed -bed12 -g $genomeSize > $sample.psl.bam
date
echo "make psl.sam"
samtools view -H $sample.psl.bam > $sample.psl.sam.header
samtools view $sample.psl.bam > $sample.psl.sam
date
echo "make psl.seq.bam"
paste $sample.psl.sam $sample.psl.query.seq.tab | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13,$11}' > $sample.psl.query.temp.seq.sam
cat $sample.psl.sam.header $sample.psl.query.temp.seq.sam > $sample.psl.query.temp.seq.header.sam
samtools view -b $sample.psl.query.temp.seq.header.sam > $sample.psl.seq.bam
date
echo "sort and index"
samtools sort -O bam -T temp.$sample $sample.psl.seq.bam > $sample.psl.seq.sort.bam
samtools index $sample.psl.seq.sort.bam
date
echo

date
echo
echo "samtools flagstat $sample.psl.seq.sort.bam"
samtools flagstat $sample.psl.seq.sort.bam
echo
echo "samtools flagstat $sample.psl.seq.circRNA.sort.bam"
samtools flagstat $sample.scan.circRNA.psl.seq.sort.bam
echo
echo
echo Done
