process VQSR_RECALIBRATE_INDEL {
    tag "VQSR-INDEL"

    input:
    path vcf
    path ref_fasta
    path ref_fai
    path ref_dict
    tuple path(mills_vcf), path(mills_tbi)
    tuple path(dbsnp_vcf), path(dbsnp_idx)

    output:
    tuple path("indel.recal"), path("indel.tranches")

    script:
    """
    gatk VariantRecalibrator \
        -R ${ref_fasta} \
        -V ${vcf} \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills_vcf} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp_idx} \
        -mode INDEL \
        --tranches-file indel.tranches \
        -O indel.recal
    """
}
