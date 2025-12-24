// atomic variant calling: one task per interval
process HAPLOTYPE_CALLER {
    tag "$sample_id - $interval"

    input:
    tuple val(sample_id), path(bam), path(bai), path(interval)
    path ref
    path ref_idx
    path ref_dict

    output:
    tuple val(sample_id), path("shard_${interval.baseName}.vcf.gz")

    script:
    """
    gatk HaplotypeCaller \
        -R ${ref} \
        -I ${bam} \
        -L ${interval} \
        -O shard_${interval.baseName}.vcf.gz
    """
}
