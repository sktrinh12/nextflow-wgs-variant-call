process MARK_DUPLICATES {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bams)

    output:
    tuple val(sample_id), path("${sample_id}_dedup.bam"), path("${sample_id}_dedup.bai")

    script:
    """
    samtools merge -@ ${task.cpus} merged.bam ${bams}

    gatk MarkDuplicatesSpark \
        -I merged.bam \
        -O ${sample_id}_dedup.bam \
        -M ${sample_id}_metrics.txt

    samtools index ${sample_id}_dedup.bam
    """
}
