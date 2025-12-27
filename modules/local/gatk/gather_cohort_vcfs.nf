// use only if genotypeGVCFs ran per interval
process GATHER_COHORT_VCFS {
    input:
    path vcfs

    output:
    tuple path("cohort.raw.vcf.gz"), path("cohort.raw.vcf.gz.tbi")

    script:
    """
    gatk GatherVcfs \
        ${vcfs.collect { "-I $it" }.join(' ')} \
        -O cohort.raw.vcf.gz

    gatk IndexFeatureFile -I cohort.raw.vcf.gz
    """
}
