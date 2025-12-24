// FINAL MERGE: Combine VCF shards
process MERGE_VCFS {
    tag "$sample_id"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path("${sample_id}_final.vcf.gz"), path("${sample_id}_final.vcf.gz.tbi")

    script:
    """
    gatk MergeVcfs \\
        ${vcfs.collect{ "-I $it" }.join(' ')} \\
        -O ${sample_id}_final.vcf.gz
    """
}
