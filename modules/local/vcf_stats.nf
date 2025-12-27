process VCF_STATS {
    tag "VCF_STATS"

    publishDir "${params.outdir}/vcf_stats", mode: 'copy'

    input:
    path vcf

    output:
    path "cohort.stats"
    path "plots/*"

    script:
    """
    bcftools stats ${vcf} > cohort.stats
    plot-vcfstats cohort.stats -p plots
    """
}
