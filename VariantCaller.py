#!/usr/bin/env python

import sys, os, shutil
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import multiprocessing as mp
import gzip
import glob
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
     __     __         _             _    ____      _ _               
     \ \   / /_ _ _ __(_) __ _ _ __ | |_ / ___|__ _| | | ___ _ __     
      \ \ / / _` | '__| |/ _` | '_ \| __| |   / _` | | |/ _ \ '__|    
       \ V / (_| | |  | | (_| | | | | |_| |__| (_| | | |  __/ |       
        \_/ \__,_|_|  |_|\__,_|_| |_|\__|\____\__,_|_|_|\___|_|       
TCCTGCTCTCCTCCTCCTATATTGAAGCCGGCGCAGGAACAGGATGAACTGTATACCCGCCCCTTTCCGG
::::::::::::::: ::::::::::::          :::::::::::::::::::::::: :::::::
TCCTGCTCTCCTCCTGCTATATTGAAGC----------ACAGGATGAACTGTATACCCGCCCTTTTCCGG
               ^            |--------|                        ^       
              SNV              INDEL                         SNV      

    VariantCaller follows the 2022 gatk & bcftools best practices.    
            Variants are also phased using WhatsHap.                  

!! Please note this script uses default settings for all gatk tools !!

#### PIPELINE ####
- Map fastq reads with bwa-mem
- Process alignments (sorting, marking duplicates)
- gatk or bcftools call SNVs and INDELS
	-m DNA      = Germline variant calling
	-m RNA      = RNAseq variant calling (identical to DNA with hard filters)
	-m SOMATIC  = Somatic variant calling (useful for tumors lines)
	-m MITO     = Mitochondrial variant calling
- WhatsHap Phasing
- bcftools/bedtools Consensus Haplotyping

#### DEPENDENCIES ####
:: All tools are expected to be in your PATH
- Python 3
- gnu-parallel
- pandas & numpy
- Pigz
- BWA
- STAR
- samtools + bcftools > 1.15
- gatk4
- gatktool
- WhatsHap
- Bedtools

#### GATK HARD FILTER SETTINGS ####
https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
INFO              SNP           INDELS
QD              <  2.0          <  2.0
QUAL            < 30.0          < 30.0
SOR             >  3.0            -   
FS              > 60.0          >200.0
MQ              < 40.0            -   
MQRankSum       <-12.5            -   
ReadPosRankSum  <- 8.0          <-20.0

### BCFTOOLS FILTER SETTINGS ###
https://www.htslib.org/workflow/filter.html
Thesholds adjusted based on local depth. 
INFO              SNP              INDELS
QUAL            <=10.0               -   
MQBZ     < -(3.5+4*DP/QUAL)   < -(5+DP/20)
RPBZ     >  (3+3*DP/QUAL)            -   
RPBZ     < -(3+3*DP/QUAL)            -   
SP          > (40+DP/2)              -   
SCBZ       > (2.5+DP/30)             -   
IDV               -                < 2.0 
IMF               -                < 0.1 
RPBZ+SCBZ         -                > 9.0
** NOTE ** SNPs are not filtered by depth (DP)

#### EXAMPLE ####
VariantCaller.py -f SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz SAMPLE_Merged.fastq.gz -r REF.fasta -c 16

#### CITE #### 
https://github.com/RhettRautsaw/Bioinformatics\n""")

################# Arguments #################

parser.add_argument("-f","--fastqs",
					type=argparse.FileType('r+'), nargs='+',
					default=None,
					help="REQUIRED: Gzipped fastq reads. (default: %(default)s)")
parser.add_argument("-r","--reference",
					type=argparse.FileType('r+'),
					default=None,
					help="REQUIRED: Unzipped reference fasta for mapping (default: %(default)s)")
parser.add_argument("-m","--mode",
					type=str,
					default='DNA',
					help="REQUIRED: \'DNA\' = Germline variant calling. \'RNA\' = RNAseq variant calling. \'SOMATIC\' = Somatic/Tumor variant calling. \'MITO\' = Mitochondrial variant calling. (default: %(default)s)")
parser.add_argument("-s","--samplename",
					type=str,
					default=None,
					help="Name of sample. Will be placed within final VCF files. (default: %(default)s)")

group1=parser.add_argument_group('Variant Calling Options')
group1.add_argument("--mpileup",
					action="store_true",
					default=False,
					help="OPTIONAL: Use bcftools mpileup/call to call SNPs rather than gatk (this is much faster than gatk HaplotypeCaller and is recommended for large files. (default: %(default)s)")
group1.add_argument("--ploidy",
					type=int,
					default=2,
					help="OPTIONAL: If dealing with a polyploid you can specify the ploidy level here. (default: %(default)s)")
group1.add_argument("--rna2genome",
					action="store_true",
					default=False,
					help="OPTIONAL: Turn on if you are mapping RNAseq data to a reference genome rather than a transcriptome. (default: %(default)s)")

group2=parser.add_argument_group('gatk Filtering Options (this section is not applicable to bcftools pipeline)')
group2.add_argument("--known_sites",
					type=argparse.FileType('r+'),
					default=None,
					help="OPTIONAL: If working with a well-annotated model organism/genome, provide known SNP sites for BaseRecalibration. If not provided, SNP calls will be bootstrapped with each iteration using the highest quality calls from the previous iteration for BaseRecalibration. Applies only to -m DNA & -m SOMATIC. (default: %(default)s)")
group2.add_argument("--cnn",
					action="store_true",
					default=False,
					help="OPTIONAL: Requires --known_sites. Turn on if INITIAL filtering should be done with GATK CNNScoreVariants. (default: %(default)s)")
group2.add_argument("-b","--bootstraps",
					type=int,
					default=1,
					help="OPTIONAL: If no known sites are available (e.g., non-model organisms), you can perform bootstrapping using the highest scoring variants from the previous iteration as known_sites for BaseRecalibration to return higher confidence variant calls. This flag specifies the number of bootstraps to perform. (default: %(default)s)")
group2.add_argument("--bootcnn",
					action="store_true",
					default=False,
					help="OPTIONAL: Turn on to filter variants using GATK CNNScoreVariants during bootstrapping using the highest scoring variants from the previous iteration as known_sites. Otherwise a hard-filter will be applied during bootstrapping. (default: %(default)s)")

group3=parser.add_argument_group('Phasing')
group3.add_argument("--nophase",
					action="store_false",
					default=True,
					help="OPTIONAL: Turn off phasing of VCF files. (default: %(default)s)")
group3.add_argument("--haplotypes",
					action="store_true",
					default=False,
					help="OPTIONAL: Turn on phased haplotype fasta creation. (default: %(default)s)")

group4=parser.add_argument_group('Multiprocessing')
group4.add_argument("-c","--cpus",
					type=int,
					default=mp.cpu_count(),
					help="OPTIONAL: Number of cpus to be used in each step (default: %(default)s)")
parser.add_argument("--version", action='version', version='VariantCaller v1.0')
args=parser.parse_args()

################# Functions #################
elog = str(len(glob.glob('VariantCaller.error*.log'))+1)
logfile = open('VariantCaller.error' + elog + '.log', 'a')

def gatkBaseRecal(tmpref, tmpbam, outpre, tmpknown):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk BaseRecalibrator :::")
	sp.call("gatk BaseRecalibrator -R " + tmpref + " -I " + tmpbam + " -O " + outpre + "_recal.table --known-sites " + tmpknown, shell=True, stdout=logfile, stderr=logfile) #stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk ApplyBQSR :::")
	sp.call("gatk ApplyBQSR -R " + tmpref + " -I " + tmpbam + " -bqsr " + outpre + "_recal.table -O " + outpre + "_recal.bam", shell=True, stdout=logfile, stderr=logfile)

## Call Variants -- gatk
def gatkCall(tmpref, tmpbam, tmpploid, outpre, outpre2, tmpcpus):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk HaplotypeCaller :::")
	sp.call("gatk HaplotypeCaller -I " + tmpbam + " -O " + outpre + "_haplotypecaller.g.vcf -R " + tmpref + " -ERC GVCF -ploidy " + tmpploid + " --native-pair-hmm-threads " + tmpcpus, shell=True, stdout=logfile, stderr=logfile)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk GenotypeGVCFs :::")
	sp.call("gatk GenotypeGVCFs -R " + tmpref + " -V " + outpre + "_haplotypecaller.g.vcf -ploidy " + tmpploid + " -O " + outpre2 + "_genotypes.vcf", shell=True, stdout=logfile, stderr=logfile)

# parallel processing chromosomes -- deprecate -- good idea but it runs out of memory easily
def gatkCallp(tmpref, tmpregions, tmpbam, tmpploid, outpre, outpre2, tmpcpus):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk HaplotypeCaller :::")
	sp.call("parallel -a " + tmpregions + " -j " + tmpcpus + " --bar --colsep \'\t\' \'gatk HaplotypeCaller -L {1} -I " + tmpbam + " -O " + outpre + "_{1}.g.vcf -R " + tmpref + " -ERC GVCF -ploidy " + tmpploid + "\'", shell=True, stdout=logfile, stderr=logfile)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk CombineGVCFs :::")
	sp.call("gatk CombineGVCFs -R " + tmpref + " -V " + outpre + "_*.g.vcf -O " + outpre2 + "_haplotypecaller.g.vcf", shell=True, stdout=logfile, stderr=logfile)
	sp.call("rm " + outpre + "_*.vcf " + outpre + "_*vcf.idx", shell=True, stdout=logfile, stderr=logfile)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk GenotypeGVCFs :::")
	sp.call("gatk GenotypeGVCFs -R " + tmpref + " -V " + outpre + "_haplotypecaller.g.vcf -ploidy " + tmpploid + " -O " + outpre2 + "_genotypes.vcf", shell=True, stdout=logfile, stderr=logfile)

## Filter Variants -- gatk 
def BootFilter(tmpvar, outpre, iter):
	# Variant Filtering -- SNPs -- Median Quality Score Cutoff to select highest scoring variants as known_sites
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Bootstrap " + iter + ": gatk SelectVariants (SNPs) :::")
	sp.call("gatk SelectVariants -V " + tmpvar + " -select-type SNP -O " + outpre + "_snp_" + iter + ".vcf", shell=True, stdout=logfile, stderr=logfile)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Bootstrap " + iter + ": Creating Upper Quantile Filter Expression :::")
	sp.call("bcftools query -f \'%QD,%QUAL,%SOR,%FS,%MQ,%MQRankSum,%ReadPosRankSum\n\' " + outpre + "_snp_" + iter + ".vcf > " + outpre + "_snp-scores_" + iter + ".csv", shell=True, stdout=logfile, stderr=logfile)
	data=pd.read_csv(outpre + "_snp-scores_" + iter + ".csv", header=None, names=["QD","QUAL","SOR","FS","MQ","MQRankSum","ReadPosRankSum"]).apply(pd.to_numeric, downcast='float', errors='coerce')
	limits = data.quantile(q=0.5, axis=0).to_frame(name='value')
	limits['sign']=['<','<','>','>','<','<','<']
	limits['filter']=limits.index
	limits['expression'] = limits['filter'] + " " + limits['sign'] + " " + limits['value'].astype(str)
	filterExp = limits['expression'].str.cat(sep=" || ")
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Bootstrap " + iter + ": \'" + filterExp + "\' :::")
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Bootstrap " + iter + ": gatk VariantFiltration (SNPs) :::")
	sp.call("gatk VariantFiltration -V " + outpre + "_snp_" + iter + ".vcf -filter \'" + filterExp + "\' --filter-name \'FILTER\' -O " + outpre + "_snp-filt_" + iter + ".vcf", shell=True, stdout=logfile, stderr=logfile)
	
	# Variant Filtering -- INDELS -- Median Quality Score Cutoff to select highest scoring variants as known_sites
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Bootstrap " + iter + ": gatk SelectVariants (INDELs) :::")
	sp.call("gatk SelectVariants -V " + tmpvar + " -select-type INDEL -O " + outpre + "_indel_" + iter + ".vcf", shell=True, stdout=logfile, stderr=logfile)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Bootstrap " + iter + ": Creating Upper Quantile Filter Expression :::")
	sp.call("bcftools query -f \'%QD,%QUAL,%FS,%ReadPosRankSum\n\' " + outpre + "_indel_" + iter + ".vcf > " + outpre + "_indel-scores_" + iter + ".csv", shell=True, stdout=logfile, stderr=logfile)
	data=pd.read_csv(outpre + "_indel-scores_" + iter + ".csv", header=None, names=["QD","QUAL","FS","ReadPosRankSum"]).apply(pd.to_numeric, downcast='float', errors='coerce')
	limits = data.quantile(q=0.5, axis=0).to_frame(name='value')
	limits['sign']=['<','<','>','<']
	limits['filter']=limits.index
	limits['expression'] = limits['filter'] + " " + limits['sign'] + " " + limits['value'].astype(str)
	filterExp = limits['expression'].str.cat(sep=" || ")
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Bootstrap " + iter + ": \'" + filterExp + "\' :::")
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Bootstrap " + iter + ": gatk VariantFiltration (INDELs) :::")
	sp.call("gatk VariantFiltration -V " + outpre + "_indel_" + iter + ".vcf -filter \'" + filterExp + "\' --filter-name \'FILTER\' -O " + outpre + "_indel-filt_" + iter + ".vcf", shell=True, stdout=logfile, stderr=logfile)
	
	# Merge and Clean
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Bootstrap " + iter + ": gatk MergeVcfs & gathering PASS variants :::")
	sp.call("gatk MergeVcfs -I " + outpre + "_snp-filt_" + iter + ".vcf -I " + outpre + "_indel-filt_" + iter + ".vcf -O " + outpre + "_combo-filt_" + iter + ".vcf", shell=True, stdout=logfile, stderr=logfile)
	sp.call("grep -E \'^#|PASS\' " + outpre + "_combo-filt_" + iter + ".vcf > " + outpre + "_combo-pass_" + iter + ".vcf", shell=True, stdout=logfile, stderr=logfile)
	#sp.call("rm " + outpre + "_snp[-_]* indels_" + iter + "* " + outpre + "*snps_" + iter + "* " + outpre + "*combined_" + iter + "*", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("gatk IndexFeatureFile -I " + outpre + "_combo-pass_" + iter + ".vcf", shell=True, stdout=logfile, stderr=logfile)

def HardFilter(tmpvar, outpre):
	# Variant Filtering -- SNPs -- Recommended Hard Filters
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk SelectVariants (SNPs) :::")
	sp.call("gatk SelectVariants -V " + tmpvar + " -select-type SNP -O " + outpre + "_snp.vcf", shell=True, stdout=logfile, stderr=logfile)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk VariantFiltration (SNPs) :::")
	sp.call("gatk VariantFiltration -V " + outpre + "_snp.vcf -filter \'QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\' --filter-name \'FILTER\' -O " + outpre + "_snp-filt.vcf", shell=True, stdout=logfile, stderr=logfile)
	
	# Variant Filtering -- INDELS -- Recommended Hard Filters
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk SelectVariants (INDELs) :::")
	sp.call("gatk SelectVariants -V " + tmpvar + " -select-type INDEL -O " + outpre + "_indel.vcf", shell=True, stdout=logfile, stderr=logfile)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk VariantFiltration (INDELs) :::")
	sp.call("gatk VariantFiltration -V " + outpre + "_indel.vcf -filter \'QD < 2.0 || QUAL < 30.0 || FS > 200.0 || ReadPosRankSum < -20.0\' --filter-name \'FILTER\' -O " + outpre + "_indel-filt.vcf", shell=True, stdout=logfile, stderr=logfile)
	
	# Merge and Clean
	sp.call("grep -E \'^#|PASS\' " + outpre + "_snp-filt.vcf > " + outpre + "_snp-pass.vcf", shell=True, stdout=logfile, stderr=logfile)
	sp.call("grep -E \'^#|PASS\' " + outpre + "_indel-filt.vcf > " + outpre + "_indel-pass.vcf", shell=True, stdout=logfile, stderr=logfile)
	sp.call("gatk MergeVcfs -I " + outpre + "_snp-pass.vcf -I " + outpre + "_indel-pass.vcf -O " + outpre + "_combo-pass.vcf", shell=True, stdout=logfile, stderr=logfile)
	sp.call("gatk IndexFeatureFile -I " + outpre + "_combo-pass.vcf", shell=True, stdout=logfile, stderr=logfile)

def BootCNN(tmpvar, tmpref, tmpknown, outpre, iter):
	# Variant Filtering -- Machine Learning Methods that require known sites
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Bootstrap " + iter + ": gatk CNNScoreVariants :::")
	sp.call("gatk CNNScoreVariants -R " + tmpref + " -V " + tmpvar + " -O " + outpre + "_cnnscore_" + iter + ".vcf", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Bootstrap " + iter + ": gatk FilterVariantTranches :::")
	sp.call("gatk FilterVariantTranches --invalidate-previous-filters -V " + outpre + "_cnnscore_" + iter + ".vcf --resource " + tmpknown + " -O " + outpre + "_combo-filt_" + iter + ".vcf", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("grep -E \'^#|PASS\' " + outpre + "_combo-filt_" + iter + ".vcf > " + outpre + "_combo-pass_" + iter + ".vcf", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

def CNNFilter(tmpvar, tmpref, tmpknown, outpre):
	# Variant Filtering -- Machine Learning Methods that require known sites
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk CNNScoreVariants :::")
	sp.call("gatk CNNScoreVariants -R " + tmpref + " -V " + tmpvar + " -O " + outpre + "_cnnscore.vcf", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk FilterVariantTranches :::")
	sp.call("gatk FilterVariantTranches --invalidate-previous-filters -V " + outpre + "_cnnscore.vcf --resource " + tmpknown + " -O " + outpre + "_combo-filt.vcf", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("grep -E \'^#|PASS\' " + outpre + "_combo-filt.vcf > " + outpre + "_combo-pass.vcf", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

## Call Variants -- bcftools
def bcfCall(tmpref, tmpregions, tmpbam, tmpploid, outpre, outpre2, tmpcpus):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ " ::: bcftools mpileup/call :::")
	sp.call("parallel -a " + tmpregions + " -j " + str(tmpcpus) + " --bar --colsep \'\\t\' \'bcftools mpileup -a FORMAT/SP -d 8000 -Ou -f " + tmpref + " -r {1} " + tmpbam + " | bcftools call -mv -Ov --ploidy " + tmpploid + " -o " + outpre + "_{1}.vcf\'", shell=True, stdout=logfile, stderr=logfile)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ " ::: bcftools concat regions :::")
	sp.call("bcftools concat " + outpre + "_*.vcf -Ov -o " + outpre2 + "_genotypes.vcf", shell=True, stdout=logfile, stderr=logfile)
	sp.call("rm " + outpre + "_*.vcf", shell=True, stdout=logfile, stderr=logfile)

## Filter Variants -- bcftools
def bcfFilter(tmpvar, outpre):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ " ::: bcftools view (filter) :::")
	snp_exp = 'QUAL <= 10 || MQBZ < -(3.5+4*DP/QUAL) || RPBZ > (3+3*DP/QUAL) || RPBZ < -(3+3*DP/QUAL) || FORMAT/SP > (40+DP/2) || SCBZ > (2.5+DP/30)'
	indel_exp = 'IDV < 2 || IMF < 0.1 || MQBZ < -(5+DP/20) || RPBZ+SCBZ > 9'
	sp.call("bcftools view -e \"(TYPE=\'SNP\' && " + snp_exp + ") || (TYPE=\'INDEL\' && " + indel_exp + ")\" " + tmpvar + " -Ov -o " + outpre + "_combo-pass.vcf", shell=True)

## Phasing -- WhatsHap
def phase(tmpref, tmpregions, tmpvcf, tmpbam, tmpploid, outpre, tmpcpus):
	if(tmpploid == 2):
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ " ::: WhatsHap Phase :::")
		sp.call("parallel -a " + tmpregions + " -j " + str(tmpcpus) + " --bar --colsep \'\\t\' " + \
		"\'bcftools view " + tmpvcf + " --regions {1} > TMP_{1}.vcf; " + \
		"samtools view " + tmpbam + " {1} -b > TMP_{1}.bam; " + \
		"samtools index TMP_{1}.bam; " + \
		"whatshap phase -o TMP_{1}.phase.vcf --indels --reference " + tmpref + " TMP_{1}.vcf TMP_{1}.bam; " + \
		"rm TMP_{1}.vcf TMP_{1}.bam*\'", shell=True, stdout=logfile, stderr=logfile)
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ " ::: bcftools concat regions :::")
		sp.call("bcftools concat TMP_*.phase.vcf -Ov -o " + outpre + "_combo-phase.vcf", shell=True, stdout=logfile, stderr=logfile)
		sp.call("rm TMP_*.phase.vcf", shell=True, stdout=logfile, stderr=logfile)
	elif(tmpploid > 2):
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ " ::: WhatsHap Polyphase :::")
		sp.call("parallel -a " + tmpregions + " -j " + str(tmpcpus) + " --bar --colsep \'\\t\' " + \
		"\'bcftools view " + tmpvcf + " --regions {1} > TMP_{1}.vcf; " + \
		"samtools view " + tmpbam + " {1} -b > TMP_{1}.bam; " + \
		"samtools index TMP_{1}.bam; " + \
		"whatshap polyphase --ploidy " + str(tmpploid) + " -o TMP_{1}.phase.vcf --indels --reference " + tmpref + " TMP_{1}.vcf TMP_{1}.bam; " + \
		"rm TMP_{1}.vcf TMP_{1}.bam*\'", shell=True, stdout=logfile, stderr=logfile)
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ " ::: bcftools concat regions :::")
		sp.call("bcftools concat TMP_*.phase.vcf -Ov -o " + outpre + "_combo-phase.vcf", shell=True, stdout=logfile, stderr=logfile)
		sp.call("rm TMP_*.phase.vcf", shell=True, stdout=logfile, stderr=logfile)
	elif(tmpploid == 1):
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Cannot phase haploid. Skipping to haplotyping ::: ")
		sp.call("cp 09_combo-pass.vcf 10_combo-phase.vcf", shell=True, stdout=logfile, stderr=logfile)

################# Setup #################

if args.reference == None:
	print("Error: no reference fasta provided. Please use (-r) to specify reference or (-h) for help.")
	quit()
if args.fastqs == None:
	print("Error: no fastq files provided. Please use (-f) to provide fastqs or (-h) for help.")
	quit()
if args.samplename == None:
	print("Error: no sample name provided. Please use (-s) to provide sample name or (-h) for help.")
	quit()

fastqs_list = [os.path.abspath(fastq.name) for fastq in args.fastqs]
fastqs = ' '.join(fastqs_list)
reference = os.path.abspath(args.reference.name)
mode = args.mode
samplename = args.samplename
ploidy = args.ploidy
if(args.known_sites!=None):
	known_sites=os.path.abspath(args.known_sites)
else:
	known_sites=args.known_sites
boots = args.bootstraps
cpus = args.cpus
cpus2 = int(cpus/4)

################# Start Variant Caller #################

print("""
     __     __         _             _    ____      _ _               
     \ \   / /_ _ _ __(_) __ _ _ __ | |_ / ___|__ _| | | ___ _ __     
      \ \ / / _` | '__| |/ _` | '_ \| __| |   / _` | | |/ _ \ '__|    
       \ V / (_| | |  | | (_| | | | | |_| |__| (_| | | |  __/ |       
        \_/ \__,_|_|  |_|\__,_|_| |_|\__|\____\__,_|_|_|\___|_|       
TCCTGCTCTCCTCCTCCTATATTGAAGCCGGCGCAGGAACAGGATGAACTGTATACCCGCCCCTTTCCGG
::::::::::::::: ::::::::::::          :::::::::::::::::::::::: :::::::
TCCTGCTCTCCTCCTGCTATATTGAAGC----------ACAGGATGAACTGTATACCCGCCCTTTTCCGG
               ^            |--------|                        ^       
              SNV              INDEL                         SNV      
""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting VariantCaller...")
print("\tReads:\t\t"+ fastqs)
print("\tReference:\t"+ reference)
print("\tMode:\t\t"+ mode)
print("\tSample Name:\t" + samplename)
if(args.mpileup):
	print("\tPipeline:\tbcftools mpileup/call")
else:
	print("\tPipeline:\tgatk HaplotypeCaller/GenotypeGVCFs")
	print("\tKnown Sites:\t" + str(known_sites))
	print("\t   CNN:\t\t" + str(args.cnn))
	print("\tBootstraps:\t" + str(boots))
	print("\t   CNN:\t\t" + str(args.bootcnn))
print("\tPloidy:\t\t" + str(ploidy))
print("\tPhase:\t\t" + str(args.nophase))
print("\t   Haplotype:\t" + str(args.haplotypes))
print("\trna2genome:\t" + str(args.rna2genome))
print("\tCPUs:\t\t"+str(cpus))

################# Pre-Processing #################

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Data Pre-Processing for Variant Calling :::")
print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Creating Reference Dictionary and Indices (if not already done) :::")
if not os.path.isfile("00_ref.bed"):
	sp.call("faidx --transform bed " + reference + " > 00_ref.bed", shell=True, stdout=logfile, stderr=logfile)
if not os.path.isfile(reference + ".fai"):
	sp.call("samtools faidx " + reference, shell=True, stdout=logfile, stderr=logfile)
if not os.path.isfile(reference.split(".fasta")[0] + ".dict"):
	sp.call("gatk CreateSequenceDictionary -R " + reference, shell=True, stdout=logfile, stderr=logfile)

# Map Data to Reference
if(args.rna2genome):
	# Map RNA to Genome (Splice Aware)
	ref_index = os.path.dirname(reference) + "/GenomeDir"
	if not os.path.isdir(ref_index):
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: STAR index :::")
		sp.call("STAR --runThreadN " + str(cpus) + " --runMode genomeGenerate --genomeDir " + ref_index + " --genomeFastaFiles " + reference, shell=True, stdout=logfile, stderr=logfile)
	if not os.path.isfile("00_aln.sam"):
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: STAR :::")
		sp.call("STAR --genomeDir " + ref_index + " --runThreadN " + str(cpus) + " --readFilesIn " + fastqs + " --outFileNamePrefix 00_aln --outSAMunmapped Within --readFilesCommand zcat --outSAMattrRGline \'@RG\\tID:readgroup1\\tSM:" + output + "\\tPL:ILLUMINA\'", shell=True, stdout=logfile, stderr=logfile)
		sp.call("mv 00_alnAligned.out.sam 00_aln.sam; rm 00_aln*.out 00_aln*.tab", shell=True, stdout=logfile, stderr=logfile)
else:
	# Map to reference genome or transcriptome
	if not os.path.isfile(reference + ".amb"):
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: bwa index :::")
		sp.call("bwa index " + reference, shell=True, stdout=logfile, stderr=logfile)
	if not os.path.isfile("00_aln.sam"):
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: bwa mem :::")
		command="bwa mem -t " + str(cpus) + " -R \'@RG\\tID:readgroup1\\tSM:" + samplename + "\\tPL:ILLUMINA\' " + reference + " " + fastqs + " > 00_aln.sam"
		#print(command)
		sp.call(command, shell=True, stdout=logfile, stderr=logfile)

# Mark Duplicates + Sort
if not os.path.isfile("01_collate.bam"):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: samtools collate :::")
	sp.call("samtools collate -@ " + str(cpus) + " 00_aln.sam 01_collate", shell=True, stdout=logfile, stderr=logfile)
if not os.path.isfile("01_fixmate.bam"):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: samtools fixmate :::")
	sp.call("samtools fixmate -@ " + str(cpus) + " -m 01_collate.bam 01_fixmate.bam", shell=True, stdout=logfile, stderr=logfile)
if not os.path.isfile("01_sort.bam"):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: samtools sort :::")
	sp.call("samtools sort -@ " + str(cpus) + " 01_fixmate.bam > 01_sort.bam", shell=True, stdout=logfile, stderr=logfile)
if not os.path.isfile("01_markdups.bam"):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: samtools markdup :::")
	sp.call("samtools markdup -@ " + str(cpus) + " -m s 01_sort.bam 01_markdups.bam", shell=True, stdout=logfile, stderr=logfile)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: samtools index :::")
	sp.call("samtools index -@ " + str(cpus) + " 01_markdups.bam", shell=True, stdout=logfile, stderr=logfile)
#if not os.path.isfile("01_markdups.bam"):
	#print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk SortSam :::")
	#sp.call("gatk SortSam -I 00_aln.sam -O 01_markdups.bam -M 01_markdups.txt", shell=True, stdout=logfile, stderr=logfile)
	#print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk MarkDuplicates :::")
	#sp.call("gatk MarkDuplicates -I 00_aln.sam -O 01_markdups.bam -M 01_markdups.txt", shell=True, stdout=logfile, stderr=logfile)
	#### SPARK FAILS...ALTERNATE GATK ABOVE.
	#print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk MarkDuplicatesSpark :::")
	#sp.call("gatk MarkDuplicatesSpark -I 00_aln.sam -O 01_markdups.bam", shell=True, stdout=logfile, stderr=logfile)

# Recalibrate Mapped Reads if Possible
if not os.path.isfile("02_recal.bam"):
	if(known_sites!=None):
		gatkBaseRecal(reference, "01_markdups.bam", "02", known_sites)
	else:
		sp.call("cp 01_markdups.bam 02_recal.bam", shell=True, stdout=logfile, stderr=logfile)
		sp.call("samtools index -@ " + str(cpus) + " 02_recal.bam", shell=True, stdout=logfile, stderr=logfile)

################# DNA #################
if(mode == "DNA"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Germline Short Variant Discovery (SNPs + Indels) :::")
	if not os.path.isfile("04_genotypes.vcf"):
		if(args.mpileup):
			bcfCall(reference, "00_ref.bed", "02_recal.bam", str(ploidy), "03","04", str(cpus))
		else:
			#gatkCallp(reference, "00_ref.bed", "02_recal.bam", str(ploidy), "03", "04", str(cpus2))
			gatkCall(reference, "02_recal.bam", str(ploidy), "03", "04", str(cpus))
	
	if not os.path.isfile("09_combo-pass.vcf") and not os.path.isfile("09_combo-pass.vcf.gz"):
		# Variant Filtering
		if(args.mpileup):
			bcfFilter("04_genotypes.vcf", "09")
			sp.call("cp 02_recal.bam 06_recal.bam", shell=True, stdout=logfile, stderr=logfile)
			sp.call("samtools index -@ " + str(cpus) + " 06_recal.bam", shell=True, stdout=logfile, stderr=logfile)
		else:
			if(known_sites!=None and args.cnn):
				CNNFilter("04_genotypes.vcf", reference, known_sites, "09", str(0))
			elif(known_sites!=None):
				HardFilter("04_genotypes.vcf", "09")
			else:
				BootFilter("04_genotypes.vcf", "05", str(0))
				# Bootstrap
				for iter in range(1, boots+1):
					prev = str(iter-1)
					iter2 = str(iter)
					print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ " ::: Bootstrap " + iter2 + " :::")
					# Base Recalibration
					gatkBaseRecal(reference, "02_recal.bam", "06", "05_combo-pass_" + prev + ".vcf")
					# Re-call Genotypes
					if(args.mpileup):
						bcfCall(reference, "00_ref.bed", "06_recal.bam", str(ploidy), "07","08", str(cpus))
					else:
						gatkCall(reference, "06_recal.bam", str(ploidy), "07", "08", str(cpus))
						#gatkCallp(reference, "00_ref.bed", "02_recal.bam", str(ploidy), "03", "04", str(cpus2))
					# Variant Filtering
					if(iter==boots): # If on the final bootstrap, apply final filter
						if(args.bootcnn):
							CNNFilter("08_genotypes.vcf", reference, "05_combo-pass_" + prev + ".vcf", "09")
						else:
							HardFilter("08_genotypes.vcf", "09")
					else: # Else apply a more stringent filtering to select only the highest scoring variants
						if(args.bootcnn):
							BootCNN("08_genotypes.vcf", reference, "05_combo-pass_" + prev + ".vcf", "05", iter2)
						else:
							BootFilter("08_genotypes.vcf", "05", iter2)

################# RNA #################
if(mode == "RNA"):
	#if(args.rna2genome):
	#	sp.call("gatk SplitNCigarReads", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: RNAseq Short Variant Discovery (SNPs + Indels) :::")
	if not os.path.isfile("04_genotypes.vcf"):
		if(args.mpileup):
			bcfCall(reference, "00_ref.bed", "02_recal.bam", str(ploidy), "03","04", str(cpus))
		else:
			#gatkCallp(reference, "00_ref.bed", "02_recal.bam", str(ploidy), "03", "04", str(cpus2))
			gatkCall(reference, "02_recal.bam", str(ploidy), "03", "04", str(cpus))
	
	if not os.path.isfile("09_combo-pass.vcf") and not os.path.isfile("09_combo-pass.vcf.gz"):
		if(args.mpileup):
			bcfFilter("04_genotypes.vcf", "09")
			sp.call("cp 02_recal.bam 06_recal.bam", shell=True, stdout=logfile, stderr=logfile)
			sp.call("samtools index -@ " + str(cpus) + " 06_recal.bam", shell=True, stdout=logfile, stderr=logfile)
		else:
			# Variant Filtering
			BootFilter("04_genotypes.vcf", "05", str(0))
			# Bootstrap
			for iter in range(1, boots+1):
				prev = str(iter-1)
				iter2 = str(iter)
				print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ " ::: Bootstrap " + iter2 + " :::")
				# Base Recalibration
				gatkBaseRecal(reference, "02_recal.bam", "06", "05_combo-pass_" + prev + ".vcf")
				# Re-call Genotypes
				#gatkCallp(reference, "00_ref.bed", "02_recal.bam", str(ploidy), "03", "04", str(cpus2))
				gatkCall(reference, "06_recal.bam", str(ploidy), "07", "08", str(cpus))
				# Variant Filtering
				if(iter==boots): # If on the final bootstrap, apply final filter
					HardFilter("08_genotypes.vcf", "09")
				else: # Else apply a more stringent filtering to select only the highest scoring variants
					BootFilter("08_genotypes.vcf", "05", iter2)

############### SOMATIC ################
if(mode == "SOMATIC"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Somatic Short Variant Discovery (SNPs + Indels) :::")
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: SORRY THIS IS NOT YET READY ::: ")
	quit()
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk Mutect2 :::")
	sp.call("gatk Mutect2", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	
	# Calculate Contamination
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk GetPileupSummaries :::")
	sp.call("gatk GetPileupSummaries", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk CalculateContamination :::")
	sp.call("gatk CalculateContamination", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	
	# Learn Orientation Bias Artifacts
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk LearnReadOrientationModel :::")
	sp.call("gatk LearnReadOrientationModel", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	
	# Filter Variants
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: gatk FilterMutectCalls :::")
	sp.call("gatk FilterMutectCalls", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

################# MITO #################
if(mode == "MITO"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Mitochondrial Short Variant Discovery (SNPs + Indels) :::")
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: SORRY THIS IS NOT YET READY ::: ")
	quit()
	sp.call("gatk PrintReads", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("gatk RevertSam", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("gatk MergeBamAlignment", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("gatk MarkDuplicates", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("gatk CollectWgsMetrics", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("gatk Mutect2 (CallMt & CallShiftedMt)", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("gatk LiftoverVcf", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("gatk MergeVcfs", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("gatk Mutect2 (MergeMutectStats)", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("gatk FilterMutectCalls (--mitochondrial-mode)", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("gatk VariantFiltration", shell=True)#, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

################# WHATSHAP #################
if(args.nophase):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: WhatsHap Phasing (SNPs + Indels) :::")
	if not os.path.isfile("10_combo-phase.vcf.gz"):
		if not os.path.isfile("09_combo-pass.vcf.gz.tbi"):
			sp.call("bgzip -@ " + str(cpus) + " 09_combo-pass.vcf", shell=True, stdout=logfile, stderr=logfile)
			sp.call("tabix 09_combo-pass.vcf.gz", shell=True, stdout=logfile, stderr=logfile)
		
		phase(reference, "00_ref.bed", "09_combo-pass.vcf.gz", "06_recal.bam", ploidy, "10", str(cpus))
		#if(ploidy == 2):
		#	sp.call("whatshap phase -o 10_combo-phase.vcf --indels --reference " + reference + " 09_combo-pass.vcf 06_recal.bam", shell=True, stdout=logfile, stderr=logfile)
		#elif(ploidy > 2):
		#	sp.call("whatshap polyphase --ploidy " + str(ploidy) + " -o 10_combo-phase.vcf --indels --reference " + reference + " 09_combo-pass.vcf 06_recal.bam", shell=True, stdout=logfile, stderr=logfile)
		#elif(ploidy == 1):
		#	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Cannot phase haploid. Skipping to haplotyping ::: ")
		#	sp.call("cp 09_combo-pass.vcf 10_combo-phase.vcf", shell=True, stdout=logfile, stderr=logfile)
		#sp.call("gatk IndexFeatureFile -I 10_combo-phase.vcf", shell=True, stdout=logfile, stderr=logfile)
		
		sp.call("bgzip -@ " + str(cpus) + " 10_combo-phase.vcf", shell=True, stdout=logfile, stderr=logfile)
		sp.call("tabix 10_combo-phase.vcf.gz", shell=True, stdout=logfile, stderr=logfile)
		sp.call("gatk SelectVariants -V 10_combo-phase.vcf.gz -select-type SNP -O 10_snp-phase.vcf", shell=True, stdout=logfile, stderr=logfile)
		sp.call("gatk SelectVariants -V 10_combo-phase.vcf.gz -select-type INDEL -O 10_indel-phase.vcf", shell=True, stdout=logfile, stderr=logfile)
	################# HAPLOTYPE FASTAS #################
	if(args.haplotypes):
		if not os.path.isfile("11_haplotype2.fasta"):
			print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Gathering Masked Haplotype Fastas :::")
			sp.call("bedtools coverage -a 00_ref.bed -b 06_recal.bam -d > 11_coverage.txt", shell=True, stdout=logfile, stderr=logfile)
			sp.call("bedtools genomecov -bga -ibam 06_recal.bam -g 00_ref.bed > 11_tmp.bed", shell=True, stdout=logfile, stderr=logfile)
			sp.call("grep -w 0$ 11_tmp.bed > 11_0cov.bed", shell=True, stdout=logfile, stderr=logfile)
			sp.call("rm 11_tmp.bed", shell=True, stdout=logfile, stderr=logfile)
			sp.call("bcftools consensus -H 1 -f " + reference + " 10_combo-phase.vcf.gz > 11_haplotmp1.fasta", shell=True, stdout=logfile, stderr=logfile)
			sp.call("bcftools consensus -H 2 -f " + reference + " 10_combo-phase.vcf.gz > 11_haplotmp2.fasta", shell=True, stdout=logfile, stderr=logfile)
			sp.call("bedtools maskfasta -fi 11_haplotmp1.fasta -bed 11_0cov.bed -fo 11_haplotype1.fasta -mc -", shell=True, stdout=logfile, stderr=logfile)
			sp.call("bedtools maskfasta -fi 11_haplotmp2.fasta -bed 11_0cov.bed -fo 11_haplotype2.fasta -mc -", shell=True, stdout=logfile, stderr=logfile)
			sp.call("perl -pi -e 's/>(.*)$/>$1\|h1/g' 11_haplo*1.fasta", shell=True, stdout=logfile, stderr=logfile)
			sp.call("perl -pi -e 's/>(.*)$/>$1\|h2/g' 11_haplo*2.fasta", shell=True, stdout=logfile, stderr=logfile)

################# CLEANUP & NAME FINAL OUTPUT FILES #################
#print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Cleaning Up :::")
#sp.call("mv 06_recal.bam " + samplename + "_reads.bam", shell=True, stdout=logfile, stderr=logfile)
#sp.call("mv 06_recal.bai " + samplename + "_reads.bai", shell=True, stdout=logfile, stderr=logfile)
#sp.call("mv 09_combo-pass.vcf " + samplename + "_unphased.vcf", shell=True, stdout=logfile, stderr=logfile)
#sp.call("mv 09_combo-pass.vcf.idx " + samplename + "_unphased.vcf.idx", shell=True, stdout=logfile, stderr=logfile)
#if(args.nophase):
#	sp.call("mv 10_combo-phase.vcf.gz " + samplename + "_phased.vcf.gz", shell=True, stdout=logfile, stderr=logfile)
#	sp.call("mv 10_combo-phase.vcf.idx " + samplename + "_phased.vcf.idx", shell=True, stdout=logfile, stderr=logfile)
#	sp.call("mv 10_combo-phase.vcf.gz.tbi " + samplename + "_phased.vcf.gz.tbi", shell=True, stdout=logfile, stderr=logfile)
#if(args.haplotypes):
#	sp.call("mv 11_haplotype1.fasta " + samplename + "_haplotype1.fasta", shell=True, stdout=logfile, stderr=logfile)
#	sp.call("mv 11_haplotype2.fasta " + samplename + "_haplotype2.fasta", shell=True, stdout=logfile, stderr=logfile)
#
#sp.call("rm 0{0..9}* 10_* 11_*", shell=True, stdout=logfile, stderr=logfile)
#sp.call("gzip *.fasta *.vcf", shell=True, stdout=logfile, stderr=logfile)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: FINISHED! :::")