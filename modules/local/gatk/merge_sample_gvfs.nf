process MERGE_SAMPLE_GVCFS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(gvcfs)

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi")

    script:
    """
    gatk MergeVcfs \
        ${gvcfs.collect { "-I $it" }.join(' ')} \
        -O ${sample_id}.g.vcf.gz

    gatk IndexFeatureFile -I ${sample_id}.g.vcf.gz
    """
}
