

process FASTP {
    label 'process'
    publishDir params.outdir, mode: 'copy', pattern: 'fastp_trimmed/*' // publish only trimmed fastq files
    publishDir params.outdir, mode: 'copy', pattern: "${params.outputname}_fastp.json" 
    publishDir params.outdir, mode: 'copy', pattern: "${params.outputname}_fastp.html"

    input:
        tuple val(sample_id), path(x)
    
    output:
        tuple val("${sample_id}") , path("fastp_trimmed/trimmed_${params.outputname}_{1,2}.fastq.gz")
        //file("${sample_id}_fastp.json")

    script:
    def qscore_cutoff = params.ontreads ? 7 : 15 //here ontreads matters

        """
        mkdir fastp_trimmed
        fastp \
        -q $qscore_cutoff \
        -i ${x[0]} -I ${x[1]} \
        -o fastp_trimmed/trimmed_${params.outputname}_1.fastq.gz -O fastp_trimmed/trimmed_${params.outputname}_2.fastq.gz \
        -j ${params.outputname}_fastp.json \
        -h ${params.outputname}_fastp.html
        """
}


process TRIMMOMATIC {
    label 'trimmomatic'
    publishDir params.outdir 


    input:
        tuple val(sample_id), path(x)
    
    output:
        tuple val("${sample_id}") , path("trimmed_paired_{1,2}.fq.gz")

    script:
        """
        trimmomatic PE ${x[0]}  ${x[1]}   \
        trimmed_paired_1.fq.gz trimmed_unpaired_1.fq.gz  \
        trimmed_paired_2.fq.gz trimmed_unpaired_2.fq.gz   \
        ILLUMINACLIP:${params.adapter}:2:30:10:2:True   \
        SLIDINGWINDOW:4:20   \
        LEADING:3   \
        TRAILING:3   \
        MINLEN:36
        
        """
}


process FASTQC {
    label 'process'
    publishDir params.outdir 


    input:
        tuple val(sample_id), path(x)
    
    output:
        path "fastqc_${sample_id}_logs"

    script:
        """
        mkdir fastqc_${sample_id}_logs
        fastqc \
        -o fastqc_${sample_id}_logs \
        -f fastq \
        -q ${x}
        """
}


process BWA {
	publishDir "${params.outdir}/MappedRead"
	label 'process'

	input:
	file ( reference )
	file ( bwa_index )
    tuple val ( sample_id ) , path ( x )

	output:
	file "${params.outputname}.sam"

    script:
	"""
	bwa mem \
    -M \
    -R '@RG\\tID:${params.rg}\\tSM:${params.samplename}\\tPL:Illumina' \
    $reference ${x[0]} ${x[1]}  > ${params.outputname}.sam

	"""
}


process BOWTIE_INDEX {
	publishDir "${params.outdir}/reference"
	label 'bowtie'

	output:
	file "index_file.1.bt2" 
	file "index_file.2.bt2"
	file "index_file.3.bt2" 
	file "index_file.4.bt2"
	file "index_file.rev.1.bt2" 
	file "index_file.rev.2.bt2"

    script:
	"""
	gunzip -dc /data/index_file.1.bt2.gz > index_file.1.bt2
	gunzip -dc /data/index_file.2.bt2.gz > index_file.2.bt2
	gunzip -dc /data/index_file.3.bt2.gz > index_file.3.bt2
	gunzip -dc /data/index_file.4.bt2.gz > index_file.4.bt2
	gunzip -dc /data/index_file.rev.1.bt2.gz > index_file.rev.1.bt2
	gunzip -dc /data/index_file.rev.2.bt2.gz > index_file.rev.2.bt2

	"""
}

process BOWTIE {
	publishDir "${params.outdir}/MappedRead"
	label 'bowtie'

	input:
	file ( index_1 )
	file ( index_2 )
    file ( index_3 )
	file ( index_4 )
    file ( index_5 )
	file ( index_6 )
    tuple val ( sample_id ) , path ( x )

	output:
	file "KM01_bowtie.sam"

    script:
	"""
    bowtie2 \
    -p 2 \
    -x index_file  \
    -1 ${x[0]} -2 ${x[1]} \
    -S KM01_bowtie.sam

	"""
}

process SAM_TO_BAM {
    publishDir params.outdir, mode: 'copy'
    label 'process'

    input:
        file ( sam_file )

    output:
       path "${params.outputname}.bam"

    script:
    """
    samtools view \
    -b \
    -@ 12 \
    ${sam_file} > ${params.outputname}.bam
    
    """    
}


process SORTING_BAM_FILE {
    publishDir params.outdir, mode: 'copy'
    label 'process'

    input:
        file ( bam_file )

    output:
       path "${params.outputname}_sorted.bam"

    script:
    """
    picard SortSam  \
    I=${bam_file}  \
    O=${params.outputname}_sorted.bam  \
    SORT_ORDER=coordinate
    
    """   
}


process MARKDUPLICATE {
    publishDir params.outdir, mode: 'copy'
    label 'process'

    input:
        file ( sorted_bam_file )   

    output:
        path "${params.outputname}_sorted_markduplicated.bam"

    script:
    """
    picard MarkDuplicates   \
    I=${sorted_bam_file}   \
    O=${params.outputname}_sorted_markduplicated.bam  \
    M=metrics.txt   \
    REMOVE_DUPLICATES=true   \
    AS=true

    """   
}


process ADD_OR_REPLACE_READGROUPS {
    publishDir params.outdir, mode: 'copy'
    label 'process'

    input:
        file ( sorted_and_markduplicated_bam_file )   

    output:
        path "${params.outputname}_ReadGroupsAdded_sorted_markduplicated.bam"

    script:
    """
    picard AddOrReplaceReadGroups   \
    I=${sorted_and_markduplicated_bam_file}  \
    O=${params.outputname}_ReadGroupsAdded_sorted_markduplicated.bam  \
    LB=SureSelectV7  \
    PL=Illumina   \
    PU=novaseq6000  \
    SM=${sorted_and_markduplicated_bam_file.baseName}  

    """   
}

process BUILDING_BAM_INDEX {
    publishDir params.outdir, mode: 'copy'
    label 'process'

    input:
        file ( sorted_markduplicated_and_readgroups_bam_file )   

    output:
        path "${params.outputname}.bai"

    script:
    """
    picard BuildBamIndex   \
    I=${sorted_markduplicated_and_readgroups_bam_file}  \
    O=${params.outputname}.bai

    """   
}

/* Process 4: Base recalibrate to detect systematic errors in base quality scores, 
            select unique alignments and index
 */
process BASE_RECALIBRATOR {
    publishDir params.outdir, mode: 'copy'
    label 'gatk'

    input:
		file ( fasta )
        file ( sorted_markduplicated_and_readgroups_bam_file )  
		file ( dbsnp )
        file ( phase1_snp )
        file ( dbsnp_index)
        file ( phase1_snp_index )
        file ( fasta_fai )
        file ( fasta3_dict )


    output:
        path "${params.outputname}_recal_BQSR.table"

    script:
    """
    gatk BaseRecalibrator   -R ${fasta}   -I ${sorted_markduplicated_and_readgroups_bam_file}  \
    --known-sites ${dbsnp}  \
    --known-sites ${phase1_snp}  \
    -O ${params.outputname}_recal_BQSR.table

    """   
}


process APPLY_BQSR {
    publishDir params.outdir, mode: 'copy'
    label 'gatk'

    input:
        file ( fasta )
        file ( sorted_markduplicated_and_readgroups_bam_file )
        file ( sorted_markduplicated_and_readgroups_recal_BQSR_table ) 
        file ( fasta_fai )
        file ( fasta_dict )


    output:
        path "${params.outputname}_recal.bam"

    script:
    """
    gatk ApplyBQSR    -R ${fasta}   -I ${sorted_markduplicated_and_readgroups_bam_file} \
    -bqsr ${sorted_markduplicated_and_readgroups_recal_BQSR_table}  \
    -O ${params.outputname}_recal.bam

    """   
}

/*
 * Process 5: Call variants with GATK HaplotypeCaller.
 *            Calls SNPs and indels simultaneously via local de-novo assembly of 
 *            haplotypes in an active region.
 *            Filter called variants with GATK VariantFiltration.    
 */
process VARIANT_CALLING {
    publishDir params.outdir, mode: 'copy'
    label 'gatk'

    input:
        file ( sorted_markduplicated_readgroups_recal_bam_file )
        file ( fasta ) 
        file ( fasta_fai )
        file ( fasta_dict )  

    output:
        path "${params.outputname}.vcf"

    script:
    """
    gatk HaplotypeCaller    -I ${sorted_markduplicated_readgroups_recal_bam_file}    \
    -O ${params.outputname}.vcf   -R ${fasta}

    """   
}

/*
VQSR stands for Variant Quality Score Recalibration is a sophisticated filtering technique applied
on the variant callset that uses machine learning to model the technical profile of variants in a training set and uses
 that to filter out probable artifacts from the callset.
 */
process VARIANTRECALIBRATOR_SNPS {
	publishDir "${params.outdir}/VariantRecalibrator"
    label 'gatk'
	
	input:
    file ( haplotypecaller_vcf )
	file ( fasta )
    file ( hapmap )
	file ( omni )
    file ( phase1_snps )
	file ( dbsnp )
	file reference_fai
	file reference_dict
	file hapmap_index
	file omni_index
	file phase1_snps_index
    file dbsnp_index

	output:
	file "recalibrate_SNP.recal"
	file "recalibrate_SNP.tranches"
    file "recalibrate_SNP.recal.idx"

	script:
	"""
	gatk VariantRecalibrator \
	-V $haplotypecaller_vcf \
 	-R $fasta \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase1_snps \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
	-an DP \
    -an QD \
	-an FS \
    -an SOR \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
	--max-gaussians 8 \
    -O recalibrate_SNP.recal \
    --tranches-file recalibrate_SNP.tranches \

	"""
}



process VQSR_APPLY_SNP {
	publishDir "${params.outdir}/VariantRecalibrator"
	label 'gatk'
	
	input:
	file haplotypecaller_vcf
	file variantrecalibrator_recal
	file variantrecalibrator_tranches
    file variantrecalibrator_recal_index

	output:
	file "recalibrated_snps_raw_indels.vcf"
	
	script:
	"""
	gatk ApplyVQSR \
	-V $haplotypecaller_vcf \
	--recal-file $variantrecalibrator_recal \
	--tranches-file $variantrecalibrator_tranches \
	-mode SNP \
	-ts-filter-level 99.0 \
	-O recalibrated_snps_raw_indels.vcf 

	"""
}

process VARIANTRECALIBRATOR_INDELS {
	publishDir "${params.outdir}/VariantRecalibrator"
	label 'gatk'
	
	input:
    file ( recalibrated_snps_raw_indels )
	file ( fasta )
    file ( golden_indel )
	file ( dbsnp )
	file ( reference_fai )
	file ( reference_dict )
	file ( golden_indel_index )
    file ( dbsnp_index )

	output:
	file "recalibrate_INDEL.recal"
	file "recalibrate_INDEL.tranches"
    file "recalibrate_INDEL.recal.idx"
	
	script:
	"""
	gatk VariantRecalibrator \
	-V $recalibrated_snps_raw_indels \
 	-R $fasta \
	--resource:mills,known=false,training=true,truth=true,prior=12.0 $golden_indel \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
	-an QD \
    -an DP \
    -an FS \
	-an SOR \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode INDEL \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
	--max-gaussians 4 \
    -O recalibrate_INDEL.recal \
    --tranches-file recalibrate_INDEL.tranches \

	"""
}

process VQSR_APPLY_INDEL {
	publishDir "${params.outdir}/VariantRecalibrator"
	label 'gatk'
	
	input:
	file recalibrated_snps_raw_indels
	file variantrecalibrator_indel_recal
	file variantrecalibrator_indel_tranches
	file variantrecalibrator_indel_recal_idx
    tuple val(sample_id), path(x)
    

	output:
	file "${params.outputname}_recalibrated_variants.vcf"
	
	script:
	"""
	gatk ApplyVQSR \
	-V $recalibrated_snps_raw_indels \
	--recal-file $variantrecalibrator_indel_recal \
	--tranches-file $variantrecalibrator_indel_tranches \
	-mode INDEL \
	-ts-filter-level 99.0 \
	-O ${params.outputname}_recalibrated_variants.vcf

	"""
}

/* Hard-filtering consists of choosing specific thresholds for one or more annotations and throwing out any variants that
 have annotation values above or below the set thresholds.
 */
process HARD_FILTERING_STEP_1 {
    publishDir params.outdir, mode: 'copy'
    label 'gatk'

    input:
        file ( vcf_file )   

    output:
        path "${params.outputname}_snp.vcf"

    script:
    """
    gatk SelectVariants  \
    -V ${vcf_file}   \
    -select-type SNP   \
    -O ${params.outputname}_snp.vcf

    """   
}


process HARD_FILTERING_STEP_2 {
    publishDir params.outdir, mode: 'copy'
    label 'gatk'

    input:
        file ( vcf_file )   

    output:
        path "${params.outputname}_indel.vcf"

    script:
    """
    gatk SelectVariants  \
    -V ${vcf_file}   \
    -select-type INDEL  \
    -O  ${params.outputname}_indel.vcf

    """   
}


process HARD_FILTERING_STEP_3 {
    publishDir params.outdir, mode: 'copy'
    label 'gatk'

    input:
        file ( vcf_snp_file )   

    output:
        path "${params.outputname}_snp_filtered.vcf"

    script:
    """
    gatk VariantFiltration -V ${vcf_snp_file} \
    -filter "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0"  \
    --filter-name "SNP_FILTER" \
    -O ${params.outputname}_snp_filtered.vcf
    
    """   
}


process HARD_FILTERING_STEP_4 {
    publishDir params.outdir, mode: 'copy'
    label 'gatk'

    input:
        file ( vcf_indel_file )   

    output:
        path "${params.outputname}_indel_filtered.vcf"
    
    script:
    """
    gatk VariantFiltration -V ${vcf_indel_file} \
    -filter "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filter-name "INDEL_FILTER" \
    -O ${params.outputname}_indel_filtered.vcf

    """   
}


process HARD_FILTERING_STEP_5 {
    publishDir params.outdir, mode: 'copy'
    label 'gatk'

    input:
        file ( vcf_snp_filtered_file )
        file ( vcf_indel_filtered_file ) 
        file ( vcf_file )  

    output:
        path "${params.outputname}_merged.vcf"

    script:
    """
    gatk MergeVcfs -I ${vcf_snp_filtered_file} -I ${vcf_indel_filtered_file } \
    -O ${params.outputname}_merged.vcf

    """   
}


process ANNOTATION {
    publishDir params.outdir, mode: 'copy'
    label 'annotation'

    input:
        file ( merged_vcf_file )
        file ( vcf_file )
  
    output:
       path "${params.outputname}_final.hg19_multianno.vcf"
       path "${params.outputname}_final.hg19_multianno.txt"
       path "${params.outputname}_final.avinput"

    script:
    """
    /opt/annovar/table_annovar.pl ${merged_vcf_file} ${params.humandb} -buildver hg19 \
    -out ${params.outputname}_final \
    --remove \
    --thread 2 \
    -protocol refGene  \
    -operation g \
    -nastring . \
    -vcfinput 
    
    """   
}

process VCF2TSV {
    publishDir params.outdir, mode: 'copy'
    label 'process'

    input:
        file ( multianno_vcf_file )
  
    output:
        path "multianno_${params.outputname}.tsv"

    script:
    """
    vcf2tsv -g  \
    -n null_string   \
    $multianno_vcf_file > multianno_${params.outputname}.tsv
   
    """   
}
