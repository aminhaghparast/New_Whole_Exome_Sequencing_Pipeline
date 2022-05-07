# Whole-Exome-Sequencing

A Nextflow pipeline for Variant Calling Analysis with NGS RNA-Seq data based on GATK best practices.

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.01.0-brightgreen.svg)](http://nextflow.io)
[![Build Status](https://github.com/CRG-CNAG/CalliNGS-NF/actions/workflows/ci.yml/badge.svg)](https://github.com/CRG-CNAG/CalliNGS-NF/actions/workflows/ci.yml)



## Quickstart 

Install Nextflow by using the following command: 

    curl -s https://get.nextflow.io | bash 
    

Clone the repository by using the following command : 

    git clone https://github.com/aminhaghparast/Whole-Exome-Sequencing-.git

Install any of Docker or Singularity for full pipeline reproducibility.

test the pipeline on a minimal dataset with a single command : 

    nextflow run main.nf -profile docker




## Components 

WES_pipeline uses the following software components and tools: 

* fastp     0.23.2
* fastqc    0.11.9
* bwa       0.7.17
* samtools  1.15
* vcflib    1.0.3
* picard    2.26.11
* GATK      4.2.6.1
* bowtie2   2.4.5
* Annovar



## Pipeline Description

The RNA sequencing (RNA-seq) data, in additional to the expression information, can be used to obtain somatic variants present in the genes of the analysed organism. The CalliNGS-NF pipeline processes RNAseq data to obtain small variants(SNVs), single polymorphisms (SNPs) and small INDELs (insertions, deletions). The pipeline is an implementation of the GATK best practices for variant calling on RNAseq and includes all major steps of the analysis, [link](http://gatkforums.broadinstitute.org/gatk/discussion/3892/the-gatk-best-practices-for-variant-calling-on-rnaseq-in-full-detail). 

In addition to the GATK best practics, the pipeline includes steps to compare obtained SNVs with known variants and to calculate allele specific counts for the overlapped SNVs.

## Input files

The CalliNGS-NF pipeline needs as the input following files:
* RNAseq reads, `*.fastq`
* Genome assembly, `*.fa`
* Known variants, `*.vcf`
* Denylisted regions of the genome, `*.bed`

The RNAseq read file names should match the following naming convention:  *sampleID{1,2}_{1,2}.extension* 

where: 
* *sampleID* is the identifier of the sample;
* the first number **1** or **2** is the replicate ID;
* the second number **1** or **2** is the read pair in the paired-end samples;
* *extension* is the read file name extension eg. `fq`, `fq.gz`, `fastq.gz`, etc. 

example: `ENCSR000COQ1_2.fastq.gz`.

## Pipeline parameters


#### `--trimming`

* The desire method for trimming reads. the available choices are "fastp" and "trimmomatic" .
* Example:
     $ nextflow run main.nf -profile docker --alignment BWA_MEM --trimming fastp


#### `--alignment`

* The desire method for alignment. the available choices are "bowtie" and "BWA-MEM" .
* Example:
    $ nextflow run main.nf -profile docker --alignment BWA_MEM --trimming fastp
 
 
#### `--reads` 
   
* Specifies the location of the paired end reads FASTQ file(s). these reads must be in a standard format as "*{1,2}*.fastq.gz"
* By default it is set to a test data located in : `$baseDir/reads/read_{1,2}.fq.gz
Example:
   $ nextflow run main.nf --reads "/home/amin/fastq/*{1,2}*.fastq.gz" --alignment BWA_MEM --trimming fastp


#### `--bedfile` 

* Target enrichment design files showing the exoms area in the used capturing kit. it must be in '.bed' format.
* By default it is set to this location: '$baseDir/data/bed_files/S31285117_Padded.bed'.

Example:

    $ nextflow run main.nf --reads "/home/amin/fastq/*{1,2}*.fastq.gz" --bed         --alignment BWA_MEM --trimming fastp


#### `--results` 
   
* Specifies the folder where the results will be stored for the user.  
* It does not matter if the folder does not exist.
* By default is set to CalliNGS-NF's folder: `results` 

Example: 

    $ nextflow run CRG-CNAG/CalliNGS-NF --results /home/user/my_results
    

    
    
## Pipeline results

For each sample the pipeline creates a folder named `sampleID` inside the directory specified by using the `--results` command line option (default: `results`).
Here is a brief description of output files created for each sample:

file | description 
---- | ----
`Recall.bam` | somatic SNVs called from the RNAseq data
`Recall.bam.bai` | comparison of the SNVs from RNAseq data with the set of known variants
`Multianno.vcf` | SNVs that are common between RNAseq calls and known variants
`Multianno.tsv` | allele counts at a positions of SNVs (only for common SNVs)
`file.json` | a histogram plot for allele frequency (only for common SNVs)
`file.html` | a histogram plot for allele frequency (only for common SNVs)


## Schematic Outline
![Image](../master/figures/workflow.png?raw=true)

## Requirements 

* [Nextflow](https://www.nextflow.io) 20.07.1 (or later)
* Java 8 or later
* [Docker](https://www.docker.com/) 1.10 (or later) or [Singularity](http://singularity.lbl.gov) engine
* [GATK](https://gatk.broadinstitute.org/) 4.1.x 

Note: CalliNGS-NF can be used without a container engine by installing in your system all the 
required software components reported in the following section. See the included 
[Dockerfile](docker/Dockerfile) for the configuration details.
 

