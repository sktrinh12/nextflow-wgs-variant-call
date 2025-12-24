nextflow.enable.dsl = 2

include { BWA_MEM } from './modules/local/bwa_mem'
include { SPLIT_FASTQ } from './modules/local/split_fastq'
include { MARK_DUPLICATES } from './modules/local/mark_dups'
include { MERGE_VCFS } from './modules/local/merge_vcfs'
include { CREATE_INTERVALS } from './modules/local/create_intervals'
include { BASE_RECALIBRATOR } from './modules/local/base_recal'
include { APPLY_BQSR } from './modules/local/apply_bqsr'
include { COLLECT_METRICS } from './modules/local/collect_metrics'
include { HAPLOTYPE_CALLER } from './modules/local/haplotype_caller'
include { SELECT_VARIANTS } from './modules/local/select_variants'

params.samplesheet = workflow.profile == 'local' ? "test_samplesheet.csv" : "samplesheet.csv"
params.ref_dir = "reference"
params.ref_name = "hg38.fa"
params.known_sites  = "reference/Homo_sapiens_assembly38.dbsnp138.vcf"

samplesheet = file("${params.base_path}/${params.samplesheet}")
ref_fasta   = file("${params.base_path}/${params.ref_dir}/${params.ref_name}")
ref_fai     = file("${params.base_path}/${params.ref_dir}/${params.ref_name}.fai")
ref_dict    = file("${params.base_path}/${params.ref_dir}/${params.ref_name.replace('.fa', '.dict')}")
known_sites = file("${params.base_path}/${params.known_sites}")
known_idx   = file("${params.base_path}/${params.known_sites}.idx")

log.info """
WGS VARIANT CALLING - N F   P I P E L I N E
============================================
input_samplesheet       : ${samplesheet}
input_ref_dir           : ${params.ref_dir}
input_known_sites       : ${known_sites}
"""

read_pairs = Channel
    .fromPath(samplesheet)
    .splitCsv(header: true)
    .map { row ->
        tuple(
            row.sample_id,
            file("${params.base_path}/${row.read1}"),
            file("${params.base_path}/${row.read2}")
        )
    }


workflow {
    // split FASTQ files into chunks for sharding
    read_chunks = SPLIT_FASTQ(read_pairs)

    // Explode chunks into individual BWA tasks
    bwa_input = read_chunks.flatMap { sample_id, files ->

        def r1 = files.findAll { it.name =~ /_1.*\.part_.*\.fastq\.gz$/ }.sort()
        def r2 = files.findAll { it.name =~ /_2.*\.part_.*\.fastq\.gz$/ }.sort()

        assert r1.size() == r2.size() : "R1/R2 chunk count mismatch for ${sample_id}"

        r1.indices.collect { i ->
            tuple(sample_id, i + 1, r1[i], r2[i])
        }
    }

    // Align each chunk (Spot instances)
    input_for_bwa = bwa_input
        .map { t ->
            tuple(t[0], t[1], t[2], t[3], file("${params.base_path}/${params.ref_dir}"), params.ref_name)
        }

    bam_shards = input_for_bwa | BWA_MEM

    // MERGE & RECALIBRATION
    // Group all shards back by sample_id for Deduplication
    dedup_ch = MARK_DUPLICATES(bam_shards.groupTuple())

    // Generate BQSR table (one task per sample)
    recal_table_ch = BASE_RECALIBRATOR(dedup_ch, ref_fasta, ref_fai, ref_dict, known_sites, known_idx)

    // VARIANT CALLING SHARDING
    // Create the 50-100 intervals for the genome
    interval_ch = CREATE_INTERVALS(ref_dict, ref_fasta, ref_fai).flatten()

    bqsr_input = dedup_ch
        .join(recal_table_ch, by: 0)
        .combine(interval_ch)

    // This creates 50 small, calibrated BAM shards
    bqsr_shards = APPLY_BQSR(bqsr_input, ref_fasta, ref_fai, ref_dict)

    vcf_shards = HAPLOTYPE_CALLER(
        bqsr_shards,
        ref_fasta, ref_fai, ref_dict
    )

    // Run metrics in the background
    COLLECT_METRICS(dedup_ch, ref_fasta)

    final_vcf = MERGE_VCFS(vcf_shards.groupTuple())
}

workflow.onComplete {
    println ( workflow.success ? \
"""
Pipeline execution summary
---------------------------
Completed at: ${workflow.complete}
Duration    : ${workflow.duration}
workDir     : ${workflow.workDir}
exit status : ${workflow.exitStatus}
""" : """
Failed: ${workflow.errorReport}
exit status : ${workflow.exitStatus}
"""
    )
}
