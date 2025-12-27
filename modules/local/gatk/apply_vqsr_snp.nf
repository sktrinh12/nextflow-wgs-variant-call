process APPLY_VQSR_SNP {
    tag "ApplyVQSR-SNP"

    input:
    path vcf
    path ref
    tuple path(recal), path(tranches)

    output:
    path "cohort.snp.vqsr.vcf.gz"

    script:
    """
    gatk ApplyVQSR \
        -R ${ref} \
        -V ${vcf} \
        --recal-file ${recal} \
        --tranches-file ${tranches} \
        --truth-sensitivity-filter-level 99.7 \
        -mode SNP \
        -O cohort.snp.vqsr.vcf.gz
    """
}
