nextflow.enable.dsl = 2

include { BWA_MEM } from './modules/local/bwa_mem'
include { MARK_DUPLICATES } from './modules/local/mark_dups'
include { CREATE_INTERVALS } from './modules/local/gatk/create_intervals'
include { BASE_RECALIBRATOR } from './modules/local/gatk/base_recal'
include { APPLY_BQSR } from './modules/local/gatk/apply_bqsr'
include { COLLECT_METRICS } from './modules/local/gatk/collect_metrics'
include { HAPLOTYPE_CALLER } from './modules/local/gatk/haplotype_caller'
include { MERGE_SAMPLE_GVCFS } from './modules/local/gatk/merge_sample_gvfs'
include { GENOMICSDB_IMPORT } from './modules/local/gatk/genomicsdb_import'
include { GENOTYPE_GVCFS } from './modules/local/gatk/genotype_gvcfs'
include { GATHER_COHORT_VCFS } from './modules/local/gatk/gather_cohort_vcfs'
include { VQSR_RECALIBRATE_SNP } from './modules/local/gatk/vqsr_recalibrate_snp'
include { VQSR_RECALIBRATE_INDEL } from './modules/local/gatk/vqsr_recalibrate_indel'
include { APPLY_VQSR_SNP } from './modules/local/gatk/apply_vqsr_snp'
include { APPLY_VQSR_INDEL } from './modules/local/gatk/apply_vqsr_indel'
include { HARD_FILTER } from './modules/local/gatk/hard_filter'
include { FUNCOTATOR } from './modules/local/gatk/funcotator'
include { VCF_STATS } from './modules/local/vcf_stats'

params.samplesheet = (
    workflow.profile.contains('local') ||
    workflow.profile.contains('k8s')
) ? 'test_samplesheet.csv' : 'samplesheet.csv'
params.ref_dir = "reference"
params.ref_name = "hg38_chr22.fa"
params.homo_sapiens = "Homo_sapiens_assembly38.dbsnp138.vcf"
params.funcotator_data_sources = "funcotator_dataSources.v1.7.20200521g"
params.hapmap = "hapmap_3.3.hg38.vcf.gz"
params.omni = "1000G_omni2.5.hg38.vcf.gz"
params.mills = "Mills_and_1000G_gold_standard.indels.vcf.gz"
params.thousandG = "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
params.dbsnp = "${params.homo_sapiens}.gz"

samplesheet = file(params.samplesheet)
hapmap      = Channel.value(tuple(
                file("${params.base_path}/${params.ref_dir}/${params.hapmap}"),
                file("${params.base_path}/${params.ref_dir}/${params.hapmap}.tbi")
            ))
omni        = Channel.value(tuple(
                file("${params.base_path}/${params.ref_dir}/${params.omni}"),
                file("${params.base_path}/${params.ref_dir}/${params.omni}.tbi")
            ))
thousandG   = Channel.value(tuple(
                file("${params.base_path}/${params.ref_dir}/${params.thousandG}"),
                file("${params.base_path}/${params.ref_dir}/${params.thousandG}.tbi")
            ))
mills       = Channel.value(tuple(
                file("${params.base_path}/${params.ref_dir}/${params.mills}"),
                file("${params.base_path}/${params.ref_dir}/${params.mills}.tbi")
            ))
dbsnp       = Channel.value(tuple(
                file("${params.base_path}/${params.ref_dir}/${params.dbsnp}"),
                file("${params.base_path}/${params.ref_dir}/${params.dbsnp}.tbi")
            ))
ref_fasta   = file("${params.base_path}/${params.ref_dir}/${params.ref_name}")
ref_fai     = file("${params.base_path}/${params.ref_dir}/${params.ref_name}.fai")
ref_dict    = file("${params.base_path}/${params.ref_dir}/${params.ref_name.replace('.fa', '.dict')}")
homo_sapiens = file("${params.base_path}/${params.ref_dir}/${params.homo_sapiens}")
homo_sapiens_idx = file("${params.base_path}/${params.ref_dir}/${params.homo_sapiens}.idx")
funcotator_sources = file("${params.base_path}/${params.ref_dir}/${params.funcotator_data_sources}")

ref_dir = file("${params.base_path}/${params.ref_dir}")
read_pairs_ch = Channel.fromPath(samplesheet)
        .splitCsv(header:true)
        .map { row ->
            tuple(
                row.sample_id,
                file("${params.base_path}/${row.read1}"),
                file("${params.base_path}/${row.read2}")
            )
        }

log.info """
WGS VARIANT CALLING - N F   P I P E L I N E
============================================
input_samplesheet       : ${samplesheet}
input_ref_dir           : ${ref_dir}
input_homo_sapiens      : ${homo_sapiens}
"""

workflow {

    // alignment
    bwa_bams = BWA_MEM(read_pairs_ch, ref_dir, params.ref_name)

    // deduplication
    dedup_bams = MARK_DUPLICATES(bwa_bams)

    // Generate BQSR table
    recal_table = BASE_RECALIBRATOR(dedup_bams, ref_fasta, ref_fai, ref_dict, homo_sapiens, homo_sapiens_idx)

    bqsr_bams = APPLY_BQSR(dedup_bams.join(recal_table), ref_fasta, ref_fai, ref_dict)

    // Run QC metrics
    COLLECT_METRICS(bqsr_bams, ref_fasta)

    // Create genome intervals (for HC + GenomicsDB)
    intervals = CREATE_INTERVALS(ref_dict, ref_fasta, ref_fai).flatten()

    // VARIANT CALLING SHARDING
    gvcf_shards = HAPLOTYPE_CALLER(
        bqsr_bams.combine(intervals),
        ref_fasta,
        ref_fai,
        ref_dict
    )

    // Merge per-sample GVCFs
    sample_gvcfs = MERGE_SAMPLE_GVCFS(gvcf_shards.groupTuple())
    gvcf_and_index = sample_gvcfs.map { sample_id, gvcf, idx -> tuple(gvcf, idx) }
    genomicsdb_input = gvcf_and_index
        .collect()           // Collect all samples into a list
        .combine(intervals)

    // Joint genotyping
    genomicsdbs = GENOMICSDB_IMPORT(
        genomicsdb_input,
        ref_fasta,
        ref_fai,
        ref_dict
    )

    cohort_vcfs = GENOTYPE_GVCFS(
        genomicsdbs,
        ref_fasta,
        ref_fai,
        ref_dict
    )

    // Gather cohort VCF
    cohort_raw_vcf  = GATHER_COHORT_VCFS(cohort_vcfs)

    // Filtering happens only after joint genotyping
    if (params.filtering == 'vqsr') {

        snp_model = VQSR_RECALIBRATE_SNP(
                        cohort_raw_vcf,
                        ref_fasta,
                        ref_fai,
                        ref_dict,
                        hapmap,
                        omni,
                        thousandG,
                        dbsnp
        )

        snp_model.view { "SNP model output: $it" }

        snp_vcf   = APPLY_VQSR_SNP(
                        cohort_raw_vcf,
                        ref_fasta,
                        snp_model
        )
        indel_model = VQSR_RECALIBRATE_INDEL(
                        snp_vcf,
                        ref_fasta,
                        ref_fai,
                        ref_dict,
                        mills,
                        dbsnp
        )
        final_vcf   = APPLY_VQSR_INDEL(
                        snp_vcf,
                        ref_fasta,
                        indel_model
        )


    } else {

        final_vcf = HARD_FILTER(
                        cohort_raw_vcf,
                        ref_fasta,
                        ref_fai,
                        ref_dict
        )

    }

    // annotation
    annotated_vcf = FUNCOTATOR(final_vcf, ref_fasta, ref_fai, ref_dict, funcotator_sources, params.ref_name)

    // final QC
    VCF_STATS(annotated_vcf)
}

workflow.onComplete {
    println ( workflow.success ? \
"""
workDir     : ${workflow.workDir}
exit status : ${workflow.exitStatus}
""" : """
Failed: ${workflow.errorReport}
exit status : ${workflow.exitStatus}
"""
    )
}
