process GENOMICSDB_IMPORT {
    tag "GenomicsDB-${interval.baseName}"

    input:
    tuple path(gvcfs), path(indices), path(interval)
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    tuple path("genomicsdb_${interval.baseName}"), path(interval)

    script:
    def variant_args = gvcfs.collect { "--variant $it" }.join(' ')
    """
    gatk GenomicsDBImport \
       --genomicsdb-workspace-path genomicsdb_${interval.baseName} \
       ${variant_args} \
       -L ${interval} \
       --reader-threads ${task.cpus}
    """
}
