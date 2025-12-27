// the default variant filtering
// <= 5 samples, non-human/no truth sets
process CNN_FILTER {
    tag "NVScore-Filter"

    input:
    tuple path(vcf), path(vcf_idx)
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    tuple path("cohort.filtered.vcf.gz"), path("cohort.filtered.vcf.gz.tbi")

    script:
    """
    gatk NVScoreVariants \
        -V ${vcf} \
        -R ${ref_fasta} \
        -O cohort.scored.vcf.gz

    gatk FilterVariantTranches \
        -V cohort.scored.vcf.gz \
        --info-key CNN_2D \
        --snp-tranche 99.95 \
        --indel-tranche 99.4 \
        -O cohort.filtered.vcf.gz

    gatk IndexFeatureFile -I cohort.filtered.vcf.gz
    """
}
