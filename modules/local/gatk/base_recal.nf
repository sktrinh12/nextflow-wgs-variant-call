process BASE_RECALIBRATOR {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref
    path ref_idx
    path ref_dict
    path known_sites
    path known_idx

    output:
    tuple val(sample_id), path("${sample_id}_recal.table")

    script:
    """
    gatk BaseRecalibrator \
        -I ${bam} \
        -R ${ref} \
        --known-sites ${known_sites} \
        -O ${sample_id}_recal.table
    """
}
