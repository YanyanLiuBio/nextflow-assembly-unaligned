process CircularPerBasePlot {
    tag "$pair_id"
    
    publishDir "${params.path}/${params.output}"

    input:
    tuple val(pair_id), path(genbank_file), path(ber_pase)

    output:
    tuple val(pair_id), path("${pair_id}_circular_coverage_plot.png")

    script:
    """
    circular_per_base_coverage_plot.py \
        --pair-id ${pair_id} \
        --genbank ${genbank_file} \
        --coverage ${ber_pase}
    """
}
