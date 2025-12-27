process FUNCOTATOR {

    input:
    tuple path(vcf), path(vcf_idx)
    path ref_fasta
    path ref_fai
    path ref_dict
    path datasources
    val ref_name

    output:
    path "cohort.annotated.vcf.gz"

    script:
    def ref_version = ref_name.replace('.fa', '').replace('_chr22', '')
    """
    gatk Funcotator \
        --variant ${vcf} \
        --reference ${ref_fasta} \
        --ref-version ${ref_version}\
        --data-sources-path ${datasources} \
        --output cohort.annotated.vcf.gz \
        --output-file-format VCF
    """
}
