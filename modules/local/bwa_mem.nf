process BWA_MEM {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2)
    path ref_dir
    val ref_name

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    bwa mem -t ${task.cpus} \\
        -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tSM:${sample_id}" \\
        ${ref_dir}/${ref_name} ${r1} ${r2} | samtools view -Sb - > ${sample_id}.bam
    """
}
