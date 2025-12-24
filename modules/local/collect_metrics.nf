process COLLECT_METRICS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta

    output:
    path "*_metrics.txt"
    path "*_histogram.pdf"

    script:
    """
    gatk CollectAlignmentSummaryMetrics \
        R=${ref} \
        I=${bam} \
        O=alignment_metrics.txt

    gatk CollectInsertSizeMetrics \
        INPUT=${bam} \
        OUTPUT=insert_size_metrics.txt \
        HISTOGRAM_FILE=insert_size_histogram.pdf
    """
}
