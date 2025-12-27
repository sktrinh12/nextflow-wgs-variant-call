process GENOTYPE_GVCFS {
    tag "${interval.baseName}"

    input:
    tuple path(genomicsdb), path(interval)
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    path "cohort.${interval.baseName}.vcf.gz"

    script:
    """
    gatk GenotypeGVCFs \
        -R ${ref_fasta} \
        -V gendb://${genomicsdb} \
        -L ${interval} \
        -O cohort.${interval.baseName}.vcf.gz
    """
}
