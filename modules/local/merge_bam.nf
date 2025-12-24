process MERGE_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned", mode: 'move'

    input:
    path sams

    output:
    path "${sample_id}.bam"

    script:
    """
    samtools merge -f ${sample_id}.bam ${sams.join(' ')}
    """
}
