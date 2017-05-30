#!/bin/bash

#############
# What: This set of bash commands downloads initial data and does some processing 
#
#		For all Individuals, chromosomes, and population:
#
# 		VCFmergeGTF
#		This code is merging the genotype calls (from vcf) with position of 
#		protein coding genes (from gencode)
#
# 		ANGSD must be verison 0.563 for the -f (filter) opition 
#		Use ANGSD to get counts, mapq filter 40, remove duplicates from BAM files 
#
# 		Getliners
#		Merge the tmp.keys (positions) with	
#		tmp.het (heterozygous protein coding SNPs for individual)
#		greps the set of keys in tmp.keys in column 2 of the file tmp.het
#
# 		ieatgor
# 		Filter for alignability output chr$chr.ind$ind.data
#		greps any entry in tmp.data that has chromosome and position within the regions 
#		specified in the "targetfile" (regions of the genome that are "callable")
#############

ind=$1
chr=$2
pop=$3
outdir=$4
rawdir=/home/user/ASEofBases/3_raw/
angsd=/home/user/ASEofBases/2_prog/angsd/angsd

echo "Results in directory:" $outdir
echo "Population: " $pop
echo "Chromosome: " $chr
echo "Individual: " $ind

### PROGRAMS
VCFmergeGTF=/home/user/ASEofBases/2_prog/VCFmergeGTF3
GETLINERS=/home/user/ASEofBases/2_prog/getliners
IEATGOR=/home/user/ASEofBases/2_prog/ieatgor

#location of BAM files
BAMSfile=$rawdir/geuvadis.bamlist

#online location of Geuvadis data
location=http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/

# merge VCF and GTF files for each individual
# note that VCFmergeGTF will only work on a phased VCF file, as it is. 
echo $ind > $outdir/tmp.ind
$VCFmergeGTF $rawdir/gencode.chr$chr.gore $rawdir/kg.chr$chr.$pop.mac2.recode.vcf $outdir/tmp.ind $outdir/tmp.merged.gz

# keep only heterozygous sites
zcat $outdir/tmp.merged.gz | awk '$8 != 0' |  awk '$8 != 3' | sed '1d' > $outdir/tmp.het
echo "Number of heterozygous protein coding SNPs for ind" $ind
wc -l $outdir/tmp.het

# list heterozygous positions
cut -f1,2 $outdir/tmp.het | sed 's/^/chr/' > $outdir/tmp.positions

### READING GEUVADIS DATA
grep $ind $BAMSfile > $outdir/tmp.bam
echo -n $location | cat - $outdir/tmp.bam > $outdir/tmp.bamlist
wc -l 4_data/tmp.bamlist

### GET COUNTS, MAPQ FILTER, REMOVES DUPLICATES SEE LAPPALAINEN SUPPLEMENT P 47
$angsd -bam $outdir/tmp.bamlist -r chr$chr: -doCounts 4 -out $outdir/tmp.counts.protein_coding.snps -dumpCounts 4 -doQsdist 1 -nThreads 15 -minMapQ 40

### COMBINING COUNTS AND POSITIONS FROM ANGSD
gunzip $outdir/tmp.counts.protein_coding.snps.pos.gz
gunzip $outdir/tmp.counts.protein_coding.snps.counts.gz
paste $outdir/tmp.counts.protein_coding.snps.pos $outdir/tmp.counts.protein_coding.snps.counts | awk '$3 != 0' | sed '1d' > $outdir/tmp.counts.protein_coding.snps.chr$chr

### FILTERING GtfVcf DATA AND MERGING WITH COUNTS
cut -f2 $outdir/tmp.counts.protein_coding.snps.chr$chr > $outdir/tmp.keys
# tmp.counts.protein_coding.snps.chr$chr is the chr, position, and counts 
# tmp.keys is a list of the positions 

$GETLINERS -c 2 -k $outdir/tmp.keys -f $outdir/tmp.het > $outdir/tmp.covered.snps.info

paste $outdir/tmp.covered.snps.info $outdir/tmp.counts.protein_coding.snps.chr$chr | cut -f1-8,11-15 | sed 's/^/chr/' >  $outdir/tmp.data

### FILTER FOR ALIGNABILITY
$IEATGOR $rawdir/chr$chr.50mer.target $outdir/tmp.data > $outdir/chr$chr.ind$ind.data

### REMOVE INTERMEDIATE FILES
rm $outdir/tmp.*
