// SHARDING: Create Interval Lists (50-100 chunks)
process CREATE_INTERVALS {
    input:
    path ref_dict
    path ref_fasta
    path ref_fai

    output: path("*.interval_list")

    script:
    """
    gatk SplitIntervals -R ${ref_fasta} -O . --scatter-count 50
    """
}
