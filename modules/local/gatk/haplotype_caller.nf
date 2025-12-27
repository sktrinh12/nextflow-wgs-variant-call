// atomic variant calling: one task per interval
process HAPLOTYPE_CALLER {
    tag "$sample_id - ${interval.baseName}"

    input:
    tuple val(sample_id), path(bam), path(bai), path(interval)
    path ref
    path ref_fai
    path ref_dict

    output:
    tuple val(sample_id), path("${sample_id}.${interval.baseName}.g.vcf.gz")

    script:
    """
    gatk HaplotypeCaller \
        -R ${ref} \
        -I ${bam} \
        -L ${interval} \
        -ERC GVCF \
        -O ${sample_id}.${interval.baseName}.g.vcf.gz
    """
}
