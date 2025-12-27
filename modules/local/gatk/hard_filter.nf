process HARD_FILTER {

    input:
    tuple path(vcf), path(vcf_idx)
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    tuple path("cohort.filtered.vcf.gz"), path("cohort.filtered.vcf.gz.tbi")

    script:
    """
    # Filter SNPs
    gatk SelectVariants \
        -R ${ref_fasta} \
        -V ${vcf} \
        --select-type-to-include SNP \
        -O snps.vcf.gz

    gatk VariantFiltration \
        -R ${ref_fasta} \
        -V snps.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "SNP_FILTER" \
        -O snps.filtered.vcf.gz

    # Filter INDELs
    gatk SelectVariants \
        -R ${ref_fasta} \
        -V ${vcf} \
        --select-type-to-include INDEL \
        -O indels.vcf.gz

    gatk VariantFiltration \
        -R ${ref_fasta} \
        -V indels.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        --filter-name "INDEL_FILTER" \
        -O indels.filtered.vcf.gz

    # Merge back
    gatk MergeVcfs \
        -I snps.filtered.vcf.gz \
        -I indels.filtered.vcf.gz \
        -O cohort.filtered.vcf.gz

    gatk IndexFeatureFile -I cohort.filtered.vcf.gz
    """
}
