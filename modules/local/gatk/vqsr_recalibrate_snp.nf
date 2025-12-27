/*
Resource	                                Purpose
hapmap.vcf.gz	                            SNP truth
omni.vcf.gz	SNP                             training
1000G.vcf.gz	                            SNP training
dbsnp.vcf.gz	                            Known variants
mills_and_1000G_gold_standard.indels.vcf.gz	INDEL truth
*/

process VQSR_RECALIBRATE_SNP {
    tag "VQSR-SNP"

    input:
    tuple path(vcf), path(vcf_idx)
    path ref_fasta
    path ref_fai
    path ref_dict

    tuple path(hapmap_vcf),   path(hapmap_tbi)
    tuple path(omni_vcf),     path(omni_tbi)
    tuple path(thousandG_vcf),path(thousandG_tbi)
    tuple path(dbsnp_vcf),    path(dbsnp_tbi)

    output:
    tuple path("snp.recal"), path("snp.tranches")

    script:
    """
    gatk VariantRecalibrator \
        -R ${ref_fasta} \
        -V ${vcf} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap_vcf} \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni_vcf} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${thousandG_vcf} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp_vcf} \
        -an QD \
        -an MQ \
        -an MQRankSum \
        -an ReadPosRankSum \
        -an FS \
        -an SOR \
        -mode SNP \
        --tranches-file snp.tranches \
        -O snp.recal
    """
}
