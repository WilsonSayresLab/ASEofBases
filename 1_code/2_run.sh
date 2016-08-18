#!/bin/bash

#############
# What: This bash script will call another script and together will conduct filtering and 
#		parsing of the genotype and RNAseq data. 
#
# 		For all Individuals, chromosomes, and population calls script 3_makeData.sh
#############

#define directories
progdir=~/ASEofBases-master/2_prog
rawdir=~/ASEofBases-master/3_raw
outdir=~/ASEofBases-master/4_data

#define input filtering files
INFOfile=$rawdir/geuvadis.E-GEUV-1.sdrf.txt
KGINDS=$rawdir/kg.inds
GTFfile=$rawdir/gencode.protein_coding.genes.v18.gtf
TRGTfile=$rawdir/wgEncodeCrgMapabilityAlign50mer.target
pop=CEU #pop= population, Central European Utah

#define program
MAKEDATA=~/ASEofBases-master/1_code/3_makeData.sh


#load modules needed for this analysis
angsd=$progdir/angsd/0.563
vcftools=$progdir/vcftools_0.1.12b/bin/vcftools

#for each autosome
for chr in {1..22}
do

	$vcftools --gzvcf $rawdir/kg.ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep $rawdir/geuvadis.$pop.kG.ind --remove-filtered-all --remove-indels --mac 2 --max-alleles 2 --recode --recode-INFO-all --out $rawdir/kg.chr$chr.$pop.mac2

# vcftools						calls the vcftools modules
# gzvcf							tells the package that the vcf file is zipped
# --keep						keep only these individuals 
# geuvadis.$pop.kG.in			list of individual sample IDs that I wish to look at, in my case the CEU population (Central European Utah) 
# --remove-filtered-all 		Removes all sites with a FILTER flag other than PASS.
# --remove-indels 				Include or exclude sites that contain an indel. For these options "indel" means any variant that alters the length of the REF allele.
# --mac 2 						Include only sites with Minor Allele Count greater than or equal to the "--mac" value and less than or equal to the "--max-mac" value. 
# --max-alleles 2 				One of these options may be used without the other. Allele count is simply the number of times that allele appears over all individuals at that site.
# --recode 						output 
# --recode-INFO-all 			These options can be used with the above recode options to define an INFO key name to keep in the output file. This option can be used multiple times to keep more of the INFO fields.	
# --out 						define output file
# kg.chr22.CEU.mac2				name of my output file 

	done
    
    #for each individual
    for ind in $indlist
    do
    $MAKEDATA $ind $chr $pop $outdir
    done

	rm $rawdir/kg.chr$chr.$pop.mac2

done

