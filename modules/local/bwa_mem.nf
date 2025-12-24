process BWA_MEM {
    tag "$sample_id - $chunk_id"

    input:
    tuple val(sample_id), val(chunk_id), path(r1), path(r2), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_${chunk_id}.bam")

    script:
    """
    bwa mem -t ${task.cpus} \\
        -R "@RG\\tID:${sample_id}_${chunk_id}\\tPL:ILLUMINA\\tSM:${sample_id}" \\
        ${ref_fasta} ${r1} ${r2} | samtools view -Sb - > ${sample_id}_${chunk_id}.bam
    """
}
