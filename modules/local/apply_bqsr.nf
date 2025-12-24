// atomic sharding: apply BQSR to small genomic intervals
process APPLY_BQSR {
    tag "$sample_id - $interval"

    input:
    tuple val(sample_id), path(bam), path(bai), path(recal_table), path(interval)
    path ref
    path ref_idx
    path ref_dict

    output:
    tuple val(sample_id), path("shard_${interval.baseName}_bqsr.bam"), path("shard_${interval.baseName}_bqsr.bai"), path(interval)

    script:
    """
    gatk ApplyBQSR \
        -I ${bam} \
        -R ${ref} \
        --bqsr-recal-file ${recal_table} \
        -L ${interval} \
        -O shard_${interval.baseName}_bqsr.bam
    """
}
