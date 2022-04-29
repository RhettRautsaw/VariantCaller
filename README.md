```
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
```

## Rhett M. Rautsaw
[![](https://img.shields.io/badge/License-GNU%20GPLv3.0-blue)](https://choosealicense.com/licenses/gpl-3.0/)

**VariantCaller** is a wrapper for the 2022 *gatk* & *bcftools* best practices + phasing with *WhatsHap*.

Generally, I recommend using the *gatk* pipeline; however, when dealing with high coverage sequencing data this pipeline can be significantly slower than the *bcftools* pipeline. Therefore, this script also facilitates *bcftools mpileup* variant calling.

!! Please note this script uses default settings or best practices for all tools !!

# Pipeline
- Map fastq reads with bwa-mem (or STAR if mapping RNA data to a genome)
- Process alignments (sorting, marking duplicates)
- gatk or bcftools call SNVs and INDELS
	- `-m DNA`		= Germline variant calling
	- `-m RNA`		= RNAseq variant calling (identical to DNA with hard filters)
	- `-m SOMATIC`	= Somatic variant calling (useful for tumors lines)
	- `-m MITO`		= Mitochondrial variant calling
- WhatsHap Phase
- bcftools/bedtools Consensus Haplotyping

# Dependencies

All of the following tools are expected to be in your PATH

- Python 3
- gnu-parallel
- pandas & numpy
- pigz
- BWA
- STAR
- samtools + bcftools > 1.15
- gatk4
- gatktool
- WhatsHap
- Bedtools

# Arguments
```
  -h, --help            show this help message and exit
  -f FASTQS [FASTQS ...], --fastqs FASTQS [FASTQS ...]
                        REQUIRED: Gzipped fastq reads. (default: None)
  -r REFERENCE, --reference REFERENCE
                        REQUIRED: Unzipped reference fasta for mapping (default: None)
  -m MODE, --mode MODE  REQUIRED: 'DNA' = Germline variant calling. 'RNA' = RNAseq variant calling. 'SOMATIC' = Somatic/Tumor variant calling. 'MITO' = Mitochondrial variant calling. (default: DNA)
  -o OUTPUT, --output OUTPUT
                        OPTIONAL: Prefix for final output files. (default: VariantCaller)
  --version             show program's version number and exit

Variant Calling Options:
  --mpileup             OPTIONAL: Use bcftools mpileup/call to call SNPs rather than gatk (this is much faster than gatk HaplotypeCaller and is recommended for large files. (default: False)
  --ploidy PLOIDY       OPTIONAL: STILL NEED TO IMPLEMENT THIS. If dealing with a polyploid you can specify the ploidy level here. (default: 2)
  --rna2genome          OPTIONAL: Turn on if you are mapping RNAseq data to a reference genome rather than a transcriptome. (default: False)

gatk Filtering Options (this section is not applicable to bcftools pipeline):
  --known_sites KNOWN_SITES
                        OPTIONAL: If working with a well-annotated model organism/genome, provide known SNP sites for BaseRecalibration. If not provided, SNP calls will be bootstrapped with each iteration using the highest quality calls from the previous iteration for BaseRecalibration. Applies only to -m DNA & -m SOMATIC. (default: None)
  --cnn                 OPTIONAL: Requires --known_sites. Turn on if INITIAL filtering should be done with GATK CNNScoreVariants. (default: False)
  -b BOOTSTRAPS, --bootstraps BOOTSTRAPS
                        OPTIONAL: If no known sites are available (e.g., non-model organisms), you can perform bootstrapping using the highest scoring variants from the previous iteration as known_sites for BaseRecalibration to return higher confidence variant calls. This flag specifies the number of bootstraps to perform. (default: 1)
  --bootcnn             OPTIONAL: Turn on to filter variants using GATK CNNScoreVariants during bootstrapping using the highest scoring variants from the previous iteration as known_sites. Otherwise a hard-filter will be applied during bootstrapping. (default: False)

Phasing:
  --nophase             OPTIONAL: Turn off phasing of VCF files. (default: True)
  --haplotypes          OPTIONAL: Turn on phased haplotype fasta creation. (default: False)

Multiprocessing:
  -c CPUS, --cpus CPUS  OPTIONAL: Number of cpus to be used in each step (default: 16)
```

# Tutorial

## gatk4
Use `gatk HaplotypeCaller/GenotypeGVCFs`
```
VariantCaller.py -f SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz SAMPLE_Merged.fastq.gz -r REF.fasta -m DNA -c 16
```

If you have known variant sites (model organism), you can supply them with `--known_sites` for gatk to recalibrate the mapped reads and reduce false positives caused by sequencing error. If no known sites are provided, then gatk will perform bootstrapping `-b` where variants are called and the highest scoring variants are used to recalibrate reads before final variants are called. 

Additionally, if you have known sites, you can use gatk's machine learning function `CNNScoreVariants` instead of hard filtering with `--cnn`. This can also be turned on in the bootstrapping step with `--bootcnn`. If these features are not specified, then the final filter will use the following settings:

Based on [gatk's 2022 best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)
| INFO			| SNPs   | INDELS |
|:-------------:|:------:|:------:|
|QD				| <  2.0 | <  2.0 |
|QUAL			| < 30.0 | < 30.0 |
|SOR			| >  3.0 |   -    |
|FS				| > 60.0 | >200.0 |
|MQ				| < 40.0 |   -    |
|MQRankSum		| <-12.5 |   -    |
|ReadPosRankSum	| <- 8.0 | <-20.0 |

## bcftools (v1.15)
Use `bcftools mpileup/call` by simply adding `--mpileup`.
```
VariantCaller.py -f SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz SAMPLE_Merged.fastq.gz -r REF.fasta -m DNA --mpileup -c 16
```

Hard filtering without bootstrapping is the only option for the `bcftools` pipeline. Filtering will use the following settings:

Based on [bcftool's recommendations](https://www.htslib.org/workflow/filter.html), thresholds are adjusted based on local depth and variant quality.
| INFO		| SNPs				| INDELS		|
|:----------|:-----------------:|:-------------:|
|QUAL		| <= 10.0			| -				|
|MQBZ		| < -(3.5+4*DP/QUAL)|  < -(5+DP/20)	|
|RPBZ		| >  (3+3*DP/QUAL)	| - 			|
|RPBZ		| < -(3+3*DP/QUAL)	| - 			|
|SP			| > (40+DP/2)		| - 			|
|SCBZ		| > (2.5+DP/30)		| - 			|
|IDV		| -					| < 2.0 		|
|IMF		| -					| < 0.1			|
|RPBZ+SCBZ	| -					| > 9.0			|

Despite recommendations, **SNPs are NOT filtered by depth (DP)**. We recommend filtering by depth at a later step with vcftools or bcftools. 

# Output

# Cite
https://github.com/RhettRautsaw/VariantCaller

Because this program only works as a wrapper for other programs, we recommend that you cite them as well. 
- [gnu-parallel]()
- [pandas & numpy]()
- [pigz]()
- [BWA]()
- [STAR]()
- [samtools + bcftools]()
- [gatk4]()
- [WhatsHap]()
- [Bedtools]()