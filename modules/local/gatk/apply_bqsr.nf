// apply BQSR to whole BAM only
process APPLY_BQSR {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai), path(recal_table)
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    tuple val(sample_id), path("${sample_id}.bqsr.bam"), path("${sample_id}.bqsr.bam.bai")

    script:
    """
    gatk ApplyBQSR \
        -I ${bam} \
        -R ${ref_fasta} \
        --bqsr-recal-file ${recal_table} \
        -O ${sample_id}.bqsr.bam

    samtools index ${sample_id}.bqsr.bam
    """
}
