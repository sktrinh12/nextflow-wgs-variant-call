process APPLY_VQSR_INDEL {
    tag "ApplyVQSR-INDEL"

    input:
    path vcf
    path ref
    tuple path(recal), path(tranches)

    output:
    path "cohort.filtered.vcf.gz"

    script:
    """
    gatk ApplyVQSR \
        -R ${ref} \
        -V ${vcf} \
        --recal-file ${recal} \
        --tranches-file ${tranches} \
        --truth-sensitivity-filter-level 99.0 \
        -mode INDEL \
        -O cohort.filtered.vcf.gz
    """
}
