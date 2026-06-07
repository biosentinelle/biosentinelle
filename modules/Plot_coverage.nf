
process Plot_coverage { // section 24.6 of the labbook 20200707
    label 'r_ig_clustering' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/figures", mode: 'copy', pattern: "{*.png}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    // publishDir "${out_path}/reports", mode: 'copy', pattern: "{plot_coverage_report.txt}", overwrite: false // 
    cache 'true'

    input:
    val fastq_name
    path cov // warning: 3 files
    path bef_read_nb
    path after_read_nb
    val ylab
    path cute_file

    output:
    path "*.png", emit: fig_ch3
    path "*.tsv", emit: cov_tsv_ch
    path "plot_coverage_report.txt", emit: plot_cov_report_ch

    script:
    """
    plot_coverage.R "${cov.baseName}" "${after_read_nb}" "${ylab}" "${fastq_name}" "${cute_file}" "plot_coverage_report.txt"
    """
}