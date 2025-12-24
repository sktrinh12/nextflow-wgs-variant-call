process BWA_MEM {
    tag "$sample_id - $chunk_id"

    input:
    tuple val(sample_id), val(chunk_id), path(r1), path(r2), path(ref_dir), val(ref_name)

    output:
    tuple val(sample_id), path("${sample_id}_${chunk_id}.bam")

    script:
    """
    bwa mem -t ${task.cpus} \\
        -R "@RG\\tID:${sample_id}_${chunk_id}\\tPL:ILLUMINA\\tSM:${sample_id}" \\
        ${ref_dir}/${ref_name} ${r1} ${r2} | samtools view -Sb - > ${sample_id}_${chunk_id}.bam
    """
}
