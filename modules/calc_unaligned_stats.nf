process CALC_UNALIGNED_STATS {
    tag "$sample_id"
    
    container 'seqwell/python:v1.0'
    input:
    tuple val(sample_id), path(unaligned_file), path(alignment_tsv)

    output:
    path("${sample_id}_unaligned_summary.csv")


    script:
    """
    calc_unaligned_stats.py ${sample_id} ${unaligned_file} ${alignment_tsv}
    """
}
