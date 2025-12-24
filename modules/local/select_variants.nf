process SELECT_VARIANTS {
    tag "select_variants"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    path vcf
    path ref

    output:
    path "raw_snps.vcf"
    path "raw_indels.vcf"

    script:
    """
    gatk SelectVariants -R ${ref} -V ${vcf} --select-type SNP   -O raw_snps.vcf
    gatk SelectVariants -R ${ref} -V ${vcf} --select-type INDEL -O raw_indels.vcf
    """
}
