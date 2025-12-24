process SPLIT_FASTQ {
    tag "$sample_id"
    publishDir "${params.outdir}/chunks", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("*.part_*.fastq.gz")

    script:
    """
    seqkit split2 -l 1000000 -1 ${read1} -2 ${read2} -O .
    """
}
