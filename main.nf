#!/usr/bin/env nextflow

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2


/* 
 * Define the default parameters
 */ 
params.outputname ="output"
params.outdir = "$baseDir/Results"
params.adapter = "$baseDir/data/adapter.fa"
params.interval = "$baseDir/data/bed_files/S31285117_Padded.bed"
params.bed = "$baseDir/data/bed_files/bed_file.bed"
params.humandb = "/opt/annovar/humandb"
params.reads = "$baseDir/reads/*{1,2}*.fastq.gz"
params.maindir = "$baseDir"
params.bam_folder = "./bam_folder"



log.info "===================================================================="
log.info "Whole Exome Sequencing Best Practice Nextflow Pipeline                        "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run oliverSI/GATK4_Best_Practice --fastq1 read_R1.fastq.gz --fastq2 read_R2.fastq.gz"
  log.info " "
  log.info "Mandatory arguments:"
  log.info "    --fastq1        FILE               Fastq(.gz) file for read1"
  log.info "    --fastq2        FILE               Fastq(.gz) file for read2"
  log.info " "
  log.info "Optional arguments:"
  log.info "    --outdir        DIR                Output directory(default: ./Results)"
  log.info "    --samplename    STRING             Sample name(dafault: fastq1 basename)"
  log.info "    --rg            STRING             Read group tag(dafault: fastq1 basename)"
  log.info " "
  log.info "===================================================================="
  exit 1
}




/* 
 * Import modules 
 */

include { 
   REFERENCE_GENOME;
   dbSNP;
   golden_indel; 
   hapmap;
   omni;
   phase1_SNPs;
   BWA_INDEX } from './modules/GATK_resource_modules.nf'


include { 
   FASTP;
   TRIMMOMATIC;
   FASTQC; 
   BOWTIE_INDEX;
   BOWTIE;
   BWA;
   SAM_TO_BAM;
   SORTING_BAM_FILE;
   MARKDUPLICATE;
   ADD_OR_REPLACE_READGROUPS;
   BUILDING_BAM_INDEX;
   BASE_RECALIBRATOR;
   APPLY_BQSR;
   VARIANT_CALLING;
   VARIANTRECALIBRATOR_SNPS;
   VQSR_APPLY_SNP;
   VARIANTRECALIBRATOR_INDELS;
   VQSR_APPLY_INDEL;
   HARD_FILTERING_STEP_1;
   HARD_FILTERING_STEP_2;
   HARD_FILTERING_STEP_3;
   HARD_FILTERING_STEP_4;
   HARD_FILTERING_STEP_5;
   ANNOTATION;
   VCF2TSV } from './modules/process_modules.nf'

/* 
 * main pipeline logic
 */

workflow {

  Channel 
      .fromFilePairs( params.reads, checkIfExists: true, size: -1 ) // default is 2, so set to -1 to allow any number of files
      .ifEmpty { error "Can not find any reads matching ${reads}" }
      .set{ read_pairs_ch }

   REFERENCE_GENOME()
  
   dbSNP()
  
   golden_indel() 
  
   hapmap()
  
   omni()
  
   phase1_SNPs()
  
   FASTQC (read_pairs_ch)


   if (!params.trimming) {exit 1, "Please specify your desirable trimming method by using --trimming argument in the command"}
   if (params.trimming== 'fastp') {
          if (!params.alignment) {exit 1, "Please specify your desirable alignment method by using --alignment argument in the command"}
          if (params.alignment== 'BWA_MEM') {
                        FASTP (read_pairs_ch) 
                        BWA_INDEX()
                        BWA (
                        REFERENCE_GENOME.out[0],
                        BWA_INDEX.out,
                        FASTP.out)
                        SAM_TO_BAM (BWA.out)
          }
          if (params.alignment== 'bowtie') {
                        FASTP (read_pairs_ch) 
                        BOWTIE_INDEX ()
                        BOWTIE (
                        BOWTIE_INDEX.out,
                        FASTP.out)
                        SAM_TO_BAM (BOWTIE.out)
          }
   }
   if (params.trimming== 'trimmomatic') {
          if (!params.alignment) {exit 1, "Please specify your desirable alignment method by using --alignment argument in the command"}
          if (params.alignment== 'BWA_MEM') {
                        TRIMMOMATIC (read_pairs_ch)
                        BWA_INDEX()
                        BWA (
                        REFERENCE_GENOME.out[0],
                        BWA_INDEX.out,
                        TRIMMOMATIC.out)
                        SAM_TO_BAM (BWA.out)
          }
          if (params.alignment== 'bowtie') {
                        TRIMMOMATIC (read_pairs_ch)
                        BOWTIE_INDEX ()
                        BOWTIE (
                        BOWTIE_INDEX.out,
                        TRIMMOMATIC.out)
                        SAM_TO_BAM (BOWTIE.out)
          }
   }

   SORTING_BAM_FILE (SAM_TO_BAM.out)

   MARKDUPLICATE (SORTING_BAM_FILE.out)

   ADD_OR_REPLACE_READGROUPS (MARKDUPLICATE.out)

   BUILDING_BAM_INDEX (ADD_OR_REPLACE_READGROUPS.out)

   BASE_RECALIBRATOR (
        REFERENCE_GENOME.out[0],
        ADD_OR_REPLACE_READGROUPS.out,
        dbSNP.out[0],
        phase1_SNPs.out[0],
        dbSNP.out[1],
        phase1_SNPs.out[1],
        REFERENCE_GENOME.out[1],
        REFERENCE_GENOME.out[2])

   APPLY_BQSR (
        REFERENCE_GENOME.out[0],
        ADD_OR_REPLACE_READGROUPS.out,
        BASE_RECALIBRATOR.out,
        REFERENCE_GENOME.out[1],
        REFERENCE_GENOME.out[2])

   VARIANT_CALLING (
        APPLY_BQSR.out,
        REFERENCE_GENOME.out[0],
        REFERENCE_GENOME.out[1],
        REFERENCE_GENOME.out[2])


   VARIANTRECALIBRATOR_SNPS (
        VARIANT_CALLING.out,
        REFERENCE_GENOME.out[0],
        hapmap.out[0],
        omni.out[0],
        phase1_SNPs.out[0],
        dbSNP.out[0],
        REFERENCE_GENOME.out[1],
        REFERENCE_GENOME.out[2],
        hapmap.out[1],
        omni.out[1],
        phase1_SNPs.out[1],
        dbSNP.out[1] )

   VQSR_APPLY_SNP (
        VARIANT_CALLING.out,
        VARIANTRECALIBRATOR_SNPS.out[0],
        VARIANTRECALIBRATOR_SNPS.out[1],
        VARIANTRECALIBRATOR_SNPS.out[2] )

   VARIANTRECALIBRATOR_INDELS (
        VQSR_APPLY_SNP.out,
        REFERENCE_GENOME.out[0],
        golden_indel.out[0],
        dbSNP.out[0],
        REFERENCE_GENOME.out[1],
        REFERENCE_GENOME.out[2],
        golden_indel.out[1],
        dbSNP.out[1] )
        
   VQSR_APPLY_INDEL (
        VQSR_APPLY_SNP.out,
        VARIANTRECALIBRATOR_INDELS.out[0],
        VARIANTRECALIBRATOR_INDELS.out[1],
        VARIANTRECALIBRATOR_INDELS.out[2],
        read_pairs_ch )

   HARD_FILTERING_STEP_1 (VARIANT_CALLING.out)

   HARD_FILTERING_STEP_2 (VARIANT_CALLING.out)

   HARD_FILTERING_STEP_3 (HARD_FILTERING_STEP_1.out)

   HARD_FILTERING_STEP_4 (HARD_FILTERING_STEP_2.out)

   HARD_FILTERING_STEP_5 (
        HARD_FILTERING_STEP_3.out,
        HARD_FILTERING_STEP_4.out,
        VARIANT_CALLING.out ) 

   ANNOTATION (
        HARD_FILTERING_STEP_5.out,
        VARIANT_CALLING.out )
        
   VCF2TSV (ANNOTATION.out[0])


}






  */

