#!/bin/bash

#############
# What: This set of bash commands downloads initial data and does some processing 
#
#               1. Convert Mapability file to bed format and create file to filter on mapability
#               2. Get protein coding gene annotations
#               3. Download & parse sample information from 1000 genomes individuals
#               4. Download & parse sample information from Geuvadis RNAseq data
#
#############

# Define directories
progdir=~/ASEofBases-master/2_prog
rawdir=~/ASEofBases-master/3_raw

#############
## 1. Convert Mapability file to bed format and create file to filter on mapability
#############
# This data is to develop a set of filters so that we can exclude regions of the 
# genome that are difficult to map to.
# 
# check to get the correct version of bigWig for your system
#
# Note: the mapability file needs to already have been downloaded an put in the BigWig directory

cd $progdir/BigWig/

# convert bigwig to bed format:
$progdir/BigWig/bigWigToBedGraph wgEncodeCrgMapabilityAlign50mer.bigWig wgEncodeCrgMapabilityAlign50mer.bed

# parse bed file to target file:
# Depending on the the system used (linux, mac, PC) you may need to change sed 's/\t/:/' to sed 's/	/:/'
#egrep "\s1$" wgEncodeCrgMapabilityAlign50mer.bed | sed 's/\t/:/' | sed 's/\t/-/' | cut -f1 > wgEncodeCrgMapabilityAlign50mer.target
egrep "\s1$" wgEncodeCrgMapabilityAlign50mer.bed | sed 's/	/:/' | sed 's/	/-/' | cut -f1 > wgEncodeCrgMapabilityAlign50mer.target
mv wgEncodeCrgMapabilityAlign50mer.target $rawdir

# parse target file to individual chromosome files
TRGTfile=$rawdir/wgEncodeCrgMapabilityAlign50mer.target
for chr in {1..22}
do
pattern="^chr$chr:"
egrep -e $pattern $TRGTfile > $rawdir/chr$chr.50mer.target
wc -l $rawdir/chr$chr.50mer.target
done

for chr in X
do
pattern="^chr$chr:"
egrep -e $pattern $TRGTfile > $rawdir/chr$chr.50mer.target
wc -l $rawdir/chr$chr.50mer.target
done

for chr in Y
do
pattern="^chr$chr:"
egrep -e $pattern $TRGTfile > $rawdir/chr$chr.50mer.target
wc -l $rawdir/chr$chr.50mer.target
done

#############
## 2. Get protein coding gene annotations
#############
# This data is to develop a set of filters so that we can analyze only variants
# in protein coding regions.
#
# gencode.v19.annotation 

cd $rawdir

# Download gencode gene annotation file
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz 

#Keep only protein coding gene entries:
zcat $rawdir/gencode.v19.annotation.gtf.gz | awk '{if($3=="gene" && $20=="\"protein_coding\";"){print $0}}' > $rawdir/gencode.protein_coding.genes.v19.gtf

### TAKING CHR FROM GTF FILE
for chr in {1..22}
do
pattern="^chr$chr"
egrep -e $pattern $rawdir/gencode.protein_coding.genes.v19.gtf | sort -n -k4 -k5 | cut -f1,4,5,9 > $rawdir/gencode.chr$chr.gore
done

for chr in X
do
pattern="^chr$chr"
egrep -e $pattern $rawdir/gencode.protein_coding.genes.v19.gtf | sort -n -k4 -k5 | cut -f1,4,5,9 > $rawdir/gencode.chr$chr.gore
done

for chr in Y
do
pattern="^chr$chr"
egrep -e $pattern $rawdir/gencode.protein_coding.genes.v19.gtf | sort -n -k4 -k5 | cut -f1,4,5,9 > $rawdir/gencode.chr$chr.gore
done

#############
## 3. Download & parse sample information from 1000 genomes individuals
#############

# Download full 1000 genomes individuals list
wget ftp://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# Make sure 1000 genomes sample list is in the raw data directory
mv integrated_call_samples_v3.20130502.ALL.panel $rawdir/kg.integrated_call_samples_v3.20130502.ALL.panel

# Cut the first column of data
cut -f1 $rawdir/kg.integrated_call_samples_v3.20130502.ALL.panel > $rawdir/kg.inds

# add the sample(s) genotype file to the 3_raw directory 
# it is not recommend to use the 1000 genomes geotype file as there may be issues with sites correctly being aligned to the Geuvidas samples 

For each chromosome, download the 1000 genomes genotype file
for chr in {1..22}
do
pattern="^chr$chr\t"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr13.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
mv ALL.chr13.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz $rawdir/kg.ALL.chr13.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
done

#############
## 4. Download & parse sample information from Geuvadis RNAseq data 
#############

# download Geuvadis sample info and move to raw data directory
wget http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt
mv E-GEUV-1.sdrf.txt $rawdir/geuvadis.E-GEUV-1.sdrf.txt

# parse to get list of bamfile names
cut -f1,29,31 $rawdir/geuvadis.E-GEUV-1.sdrf.txt | grep _1.fastq.gz | cut -f3 > $rawdir/geuvadis.bamlist

# list of individual ids
cut -f1,29,31 $rawdir/geuvadis.E-GEUV-1.sdrf.txt | grep _1.fastq.gz | cut -f1 > $rawdir/geuvadis.ind

## get the list of individual ids in each population 
## all Geuvadis samples are present in 1000 genomes
for pop in CEU FIN GBR TSI YRI
do
cut -f1,29,31,33 $rawdir/geuvadis.E-GEUV-1.sdrf.txt | grep _1.fastq.gz | grep $pop | cut -f1 > $rawdir/geuvadis.$pop.ind
grep -f $rawdir/geuvadis.$pop.ind $rawdir/kg.inds > $rawdir/geuvadis.$pop.kG.ind
done

