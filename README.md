# ASEofBases: Allele Specific Expression method and pipeline for analysis

#### What is ASEofBases?
Biased allele expression, refers to the imbalanced expression of the two alleles in a diploid genome. Unequal transcription of alleles may occur due to cis-regulatory element variation or allele-specific epigenetic modifications.  Allelic imbalance may be incorrectly inferred due to technical variation inherent in RNA-Seq data, including read depth, reference mapping bias, and the overdispersion of reads. To correct for technical variation we develop a logistic regression model with a mixed effects approach to combine information regarding biased allele expression from many individuals in a population, and across multiple genes. Here we describe a new method for inferring allelic imbalance that combines information from multiple SNP sites within a transcribed unit, using a logistic regression model that explicitly models the effects of reference bias and SNP type biases. Additionally, by using a mixed effects approach the method also makes it possible to combine information from many individuals in a population for each gene and to test hypothesis regarding allelic imbalance in specific genes within and across populations. 

#### Citing ASEofBases
If you use ASEofBases package, please cite:

- Inference on imbalanced allelic expression from RNA sequencing data in individuals and across groups using logistic regression models. Skotte L, Olney K, Nielsen R, and Wilson Sayres M. (https://dx.doi.org/10.6084/m9.figshare.1564499.v1) Manuscript in preparation

If you use the `1000Genomes` or `Geuvidas` data sets, please also cite:

-  A global reference for human genetic variation, The 1000 Genomes Project Consortium, Nature 526, 68-74 (01 October 2015) doi: [10.1038/nature15393](http://www.1000genomes.org/home).
- Lappalainen et al. Nature 2013 : Transcriptome and genome sequencing uncovers functional variation in humans (http://dx.doi.org/10.1038/nature12531).

#### Key features
- method model of allele-specific expression in single individuals  
- method model of allele-specific expression using "population data"

Additionally code includes simulations:
- simulations for modeling statistcal power for various read depth and various number of SNPs per gene 
- simulations for modeling SNP-wise binomial, binomial-based logistic regression, and binomial-based logistic regression correcting for the over dispersion of read counts. 


#### Installation
Download [the latest release of ASEofBases](https://github.com/WilsonSayresLab/ASEofBases/archive/master.zip), unzip and add the `bin` directory to your `PATH`. E.g.:

		wget -O ASEofBases-master.zip https://github.com/WilsonSayresLab/ASEofBases/archive/master.zip
		unzip ASEofBases-master.zip && cd ASEofBases-master
		export PATH=`pwd`/bin:$PATH

You will also need the following perviosuly published programs installed:

- angsd0.563				https://github.com/ANGSD/angsd
- vcftools_0.1.12b			http://sourceforge.net/projects/vcftools/files/
- zlib-1.2.8				http://www.zlib.net	
- samtools/htslib			https://github.com/samtools/htslib
- bigwig				https://genome.ucsc.edu/goldenPath/help/bigWig.html
- bcftools				https://samtools.github.io/bcftools/
- bwa 					http://bio-bwa.sourceforge.net/bwa.shtml
- GATK GenomeAnalysisTK.jar 		https://software.broadinstitute.org/gatk/documentation/quickstart
- snpEff				http://snpeff.sourceforge.net/download.html

#### Directory overview of ASEofBases:

/ASEofBases/1_code/	# bash scripts
- 1_get.sh
- 2_run.sh
- 3_makeData.sh
- ASEofBasesAnalysis.R
- simAoB.R

/ASEofBases/2_prog/	# locally written C++ programs
- getliners
- ieatgor
- VCFmergeGTF3
- perviosuly published programs which will be installed. See Step 2 program check for more information. 
 
/ASEofBases/3_raw/	# raw and parsed data
- target mapability files
- annotation file
- raw vcf file
- filtered vcf file
- individual ID list

/ASEofBases/4_data/	# parsed RNAseq data for filtering & analysis in R
- idividual ouptut data (sampleID.data)

/ASEofBases/5_out/	# analysis and output made in R
- tab delimited output files

/ASEofBases/6_plot/	# plots made from the output data made in R 
- plots in pdf format

# ASEofBases pipeline for human chromosome 13 sample data set
### Step 1. Compile programs/code to be executable
Getliners - Merge the tmp.keys (positions) with tmp.het (heterozygous protein coding SNPs for individual) and greps the set of keys in tmp.keys in column 2 of the file tmp.het

ieatgor - Filter for alignability output chr$chr.ind$ind.data, greps any entry in tmp.data that has chromosome and position within the regions and specified in the "targetfile" (regions of the genome that are "callable")

VCFmergeGTF - This code is merging the genotype calls (from vcf) with position of protein coding genes (from gencode)

		cd 2_prog/
		g++ -O3 -o VCFmergeGTF3 VCFmergeGTF3.cpp -lz
		g++ -O3 -o ieatgor ieatgorV2.cpp -lz
		g++ -O3 -o getliners getliners.cpp -lz

### Step 2. Program check
Make sure each of these programs are installed and added to the /ASEofBases/2_prog/ directory 
- bigwig (https://genome.ucsc.edu/goldenPath/help/bigWig.html)
- angsd (https://github.com/ANGSD/angsd)
- vcftools_0.1.12b (http://sourceforge.net/projects/vcftools/files/)
- zlib-1.2.8 (http://www.zlib.net)
- samtools/htslib (https://github.com/samtools/htslib)


		# angsd and samtools
		Using a local folder containing htslib
		git clone https://github.com/samtools/htslib.git;
		git clone https://github.com/angsd/angsd.git;
		cd htslib;make;cd ../angsd;make HTSSRC=../htslib
		cd ..
		
		Systemwide installation of htslib
		git clone https://github.com/angsd/angsd.git;
		cd angsd;make
		cd ..
		
		# Bcftools by samtools
		git clone --branch=develop git://github.com/samtools/htslib.git
		git clone --branch=develop git://github.com/samtools/bcftools.git
		git clone --branch=develop git://github.com/samtools/samtools.git
		cd bcftools; make
		cd ../samtools; make
		cd ..
		
		# bwa 
		wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2
		tar -vxjf bwa-0.7.15.tar.bz2
		cd bwa-0.7.15
		make
		cd ..
		
		# Picard tools 
		wget https://github.com/broadinstitute/picard/archive/master.zip
		unzip master.zip 
		cd picard-master/
		./gradlew jar
		./gradlew clean
		cd ..
		
		# snpEff
		wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
		unzip snpEff_latest_core.zip
		
		# vcftoolss_0.1.12b
		wget https://sourceforge.net/projects/vcftools/files/vcftools_0.1.12b.tar.gz
		tar -xzvf vcftools_0.1.12b.tar.gz
		cd  vcftools_0.1.12b
		make
		cd ..
		rm vcftools_0.1.12b.tar.gz
		
		# zlib
		wget http://zlib.net/zlib-1.2.8.tar.gz
		tar -xzvf zlib-1.2.8.tar.gz
		rm zlib-1.2.8.tar.gz
		
		#BigWig
		mkdir BigWig
		cd BigWig
		rsync -aP rsync://hgdownload.cse.ucsc.edu/genome/admin/exe/linux.x86_64.v287/ ./


### 3. Download Mapability file for filtering for uniqueness of the reference genome from ENCODE 
The mapability file is used for identifing uniqueness of the reference GRCh37/hg19 genome assembly. Mapability files were generated using different window sizes, high signal will be found in areas where the sequence is unique.

Download the Mapability file:
	
		wget http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig
		mv wgEncodeCrgMapabilityAlign50mer.bigWig /ASEofBases/2_prog/BigWig/

### 4. Get raw data and conduct initial parsing
1_getRaw.sh is a bash script for downloading initial data and does some processing. 
- 1_getRaw.sh Overview:
	1. Convert Mapability file to bed format and create file to filter on mapability
	2. Get protein coding gene annotations
	3. Download & parse individual genotype information from `1000Genomes` individuals
	4. Download & parse individual information RNAseq data from `Geuvadis`

			sh /ASEofBases/bash_scripts/1_getRaw.sh


- Outputs:
	1. **chr#.50mer.target** - target file to individual chromosome files. Example: chr1:10207-10212
	2. **gencode.chr#.gore** - protein coding gene entries. Example: chr1	69091	70008	gene_id "ENSG00000186092.4"; 
	3. **gencode.protien_coding.genes.v19.gtf** - gencode gene annotation file, containing only protein coding regions 
	4. **chr#.genotypes.vcf** - 1000Genomes genotype file
	5. **individual_ID_list** -  Geuvadis RNAseq individual id list

## 5. Parse and merge genotype and transcriptome data
This bash script will call another script and together will conduct filtering and parsing of the genotype and RNAseq data.

		sh /ASEofBases/1_code/2_run.sh # this program runs 3_makeData.sh

- Output: will result in a file for each individual with these columns:
	1. chr: chromosome 
	2. pos: position
	3. ref: reference allele
	4. alt: alternative allele
	5. g: gene
	6. a1: observed phased allele 1
	7. a2: observed phased allele 2
	8. type: snp type (1 if allele 1 is the alternative allele and 2 if allele 2 is the alternative allele)
	9. n: total counts of reads overlapping
	10. nA: counts of As
	11. nC: counts of Cs
	12. nG: counts of Gs
	13. nT: counts of Ts

## 6. Infer allele specific expression using ASEofBases on the filtered data that was generated in the pervious steps
For regression analysis using data generated in pervious steps

		ASEofBases.R

# ASEofBases pipeline for simulated data 
For simulations

		simAoB.R

![simulations collage](https://raw.github.com/WilsonSayresLab/ASEofBases/master/doc/simulations.png)	
