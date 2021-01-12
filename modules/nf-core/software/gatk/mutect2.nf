include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
def options    = initOptions(params.options)

environment = params.enable_conda ? "bioconda::gatk4-spark=4.1.8.1" : null
container = "quay.io/biocontainers/gatk4-spark:4.1.8.1--0"
if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) container = "https://depot.galaxyproject.org/singularity/gatk4-spark:4.1.8.1--0"

process GATK_MUTECT2 {
    label 'CPUS_2'

    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda environment
    container container

    input:
        tuple val(meta), path(bam), path(bai), file(interval)
        path germline_resource
        path germline_resource_index
        path dict
        path fasta
        path fai
        path pon
        path ponIndex

    output:
        tuple val(meta), path(interval), path("${interval.baseName}_${meta.id}.mutect2.vcf"), emit: interval_mutect2_vcf
        tuple val(meta), path("*.stats")                                                      emit: interval_stats
        tuple val(meta), path("*.stats"), path("*.vcf")                                       emit: mutect2Stats optional true

    script:
    interval_option   = params.no_intervals ? "" : "-L ${interval}"
    output_option     = params.no_intervals ? "-O ${meta.id}.mutect2.vcf" : "-O ${interval.baseName}_${meta.id}.mutect2.vcf"
    soft_clipped_option = params.ignore_soft_clipped_bases ? "--dont-use-soft-clipped-bases true" : ""
    pon = params.pon ? "--panel-of-normals ${pon}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        Mutect2 \
        -R ${fasta} \
        -I ${normal_bam} \
        -I ${tumour_bam} \
        ${interval_option} \
        ${soft_clipped_option} \
        ${output_option} \
        --germline-resource ${germline_resource} \
        ${pon}
    """
}
