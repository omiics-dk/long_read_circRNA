#!/bin/bash


# Make folders
mkdir data
mkdir data/human
mkdir data/mouse
mkdir test_fastq

# Download data
wget -O test_fastq/human_brain_100k.fq.gz https://www.dropbox.com/s/6i7gkb3r81u2kpn/human_brain_100k.fq.gz?dl=1

wget -O data/human/circAtlas2.0_June2019_human_hg19_circRNA.0-based.bed https://www.dropbox.com/s/7dmycsjn23lcly9/circAtlas2.0_June2019_human_hg19_circRNA.0-based.bed?dl=1
wget -O data/human/CIRCpedia_v2_June2019_human_hg19_All_circRNA.unique.bed https://www.dropbox.com/s/qj2h3y9s67efnlu/CIRCpedia_v2_June2019_human_hg19_All_circRNA.unique.bed?dl=1
wget -O data/human/hg19.chrom.sizes https://www.dropbox.com/s/xiag413fshhzsam/hg19.chrom.sizes?dl=1
wget -O data/human/hg19.fa https://www.dropbox.com/s/xms15iykznabq7p/hg19.fa?dl=1
wget -O data/human/hg19.fa.fai https://www.dropbox.com/s/4g2k5l1p6il2jsh/hg19.fa.fai?dl=1
wget -O data/human/hsa_circRNA_complete.hg19.unique.sort.length.bed https://www.dropbox.com/s/3z03wfnntnxvlyz/hsa_circRNA_complete.hg19.unique.sort.length.bed?dl=1
wget -O data/human/Human_refFlat_exon_hg19_Oct2018.merge.bed https://www.dropbox.com/s/125zie5qd4fg96c/Human_refFlat_exon_hg19_Oct2018.merge.bed?dl=1
wget -O data/human/Human_refFlat_exon_hg19_Oct2018.sort.bed https://www.dropbox.com/s/xt0r1xvs5r1x7fk/Human_refFlat_exon_hg19_Oct2018.sort.bed?dl=1
wget -O data/human/Human_refFlat_hg19_Oct2018.unique.merge.bed https://www.dropbox.com/s/zc0v9tb1j0v6206/Human_refFlat_hg19_Oct2018.unique.merge.bed?dl=1
wget -O data/human/UCSC-EST-exons_hg19_09-2018.bed https://www.dropbox.com/s/tktvsayp5ahyut0/UCSC-EST-exons_hg19_09-2018.bed?dl=1

wget -O data/mouse/circAtlas2.0_Aug2019_mouse_mm10_circRNA.0-based.bed https://www.dropbox.com/s/7qe6m97jcauz4cy/circAtlas2.0_Aug2019_mouse_mm10_circRNA.0-based.bed?dl=1
wget -O data/mouse/CIRCpedia_v2_June2019_mouse_mm10_All_circRNA.unique.bed https://www.dropbox.com/s/97uej0pp37j9cz1/CIRCpedia_v2_June2019_mouse_mm10_All_circRNA.unique.bed?dl=1
wget -O data/mouse/mm10.chrom.sizes https://www.dropbox.com/s/kvv21ewlxw5wbjg/mm10.chrom.sizes?dl=1
wget -O data/mouse/mmu_circRNA_complete.mm10lift.sort.length.bed https://www.dropbox.com/s/6em2zqmbq1qq9hr/mmu_circRNA_complete.mm10lift.sort.length.bed?dl=1
wget -O data/mouse/Mouse_refFlat_exon_mm10_Oct2018.merge.bed https://www.dropbox.com/s/r9eyq7b79fdywmq/Mouse_refFlat_exon_mm10_Oct2018.merge.bed?dl=1
wget -O data/mouse/Mouse_refFlat_exon_mm10_Oct2018.sort.bed https://www.dropbox.com/s/wup35acnp5p70im/Mouse_refFlat_exon_mm10_Oct2018.sort.bed?dl=1
wget -O data/mouse/Mouse_refFlat_mm10_Oct2018.unique.merge.bed https://www.dropbox.com/s/lmsrt9yo9xrd5tt/Mouse_refFlat_mm10_Oct2018.unique.merge.bed?dl=1
wget -O data/mouse/Mus_musculus.GRCm38.87.chr-fix.fa https://www.dropbox.com/s/8sojnmhc5z9kizs/Mus_musculus.GRCm38.87.chr-fix.fa?dl=1
wget -O data/mouse/UCSC-EST-exons_mm10_09-2018.bed https://www.dropbox.com/s/2lbll2mjipyp96k/UCSC-EST-exons_mm10_09-2018.bed?dl=1


# Unpack needed tools
cd scripts
tar -zxvf samtools-0.1.18.tar.gz
tar -zxvf pblat.tar.gz

echo
echo "Done!"
