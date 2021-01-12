/*
================================================================================
                            SOMATIC VARIANT CALLING
================================================================================
*/

params.mutect2_options                = [:]
params.platypus_options               = [:]


include { GATK_MUTECT2 as MUTECT2 }                 from '../../nf-core/software/gatk/mutect2' addParams(options: params.mutect2_options)
include { GATK_GENOTYPEGVCF as GENOTYPEGVCF }       from '../../nf-core/software/gatk/genotypegvcf'    addParams(options: params.genotypegvcf_options)
include { CONCAT_VCF as CONCAT_GVCF }               from '../process/concat_vcf'                       addParams(options: params.concat_gvcf_options)
include { CONCAT_VCF as CONCAT_HAPLOTYPECALLER }    from '../process/concat_vcf'                       addParams(options: params.concat_haplotypecaller_options)
include { STRELKA_GERMLINE as STRELKA }             from '../../nf-core/software/strelka/germline'     addParams(options: params.strelka_options)

workflow SOMATIC_VARIANT_CALLING {
    take:
        bam        // channel: [mandatory] bam
        dbsnp      // channel: [mandatory] dbsnp
        dbsnp_tbi  // channel: [mandatory] dbsnp_tbi
        dict       // channel: [mandatory] dict
        fai        // channel: [mandatory] fai
        fasta      // channel: [mandatory] fasta
        intervals  // channel: [mandatory] intervals
        target_bed // channel: [optional]  target_bed
        tools      //   list:  [mandatory] list of tools

    main:

    if ('mutect2' in tools) {

// separate BAM by status
        bam_normal = Channel.create()
        bam_tumor  = Channel.create()

        bam
         .choice(bam_tumor, bam_normal) {statusMap[it[0], it[1]] == 0 ? 1 : 0}
        // Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
        // Remapping channel to remove common key idPatient
        pair_bam = bam_normal.cross(bam_tumor).map {
            normal, tumor ->
            [normal[0], normal[1], tumor[1]]
            }

        pair_bam = pair_bam.dump(tag:'BAM Somatic Pair')
        mutect2_interval_bam = pair_bam.combine(intervals)
    
    if ('mutect2' in tools) {

        MUTECT2(
            mutect2_interval_bam,
            pon,
            pon_index,
            germline_resource,
            germline_resource_index,
            dict,
            fasta,
            fai,
            target_bed)

    }
        if ('platypus' in tools) {

        PLATYPUS(
            MUTECT2.out.interval_mutect2_vcf
            mutect2_interval_bam,
            pon,
            pon_index,
            germline_resource,
            germline_resource_index,
            dict,
            fasta,
            fai,
            target_bed)

        }

}
