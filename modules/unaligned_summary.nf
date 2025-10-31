process UNALIGNED_SUMMARY {
    tag "combine summaries"

    input:
    path summary_files

    output:
    path "all_samples_unaligned_summary.csv"

    script:
    """
    # Take the first header line from the first file
    head -n 1 \$(ls -1 *.csv | head -n 1) > all_samples_unaligned_summary.csv

    # Append all data lines (skip header lines)
    for f in \$(ls *contigs*.csv); do
        tail -n +2 "\$f" >> all_samples_unaligned_summary.csv
    done
    """
}