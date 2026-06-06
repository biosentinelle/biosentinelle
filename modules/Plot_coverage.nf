
process Plot_coverage { // section 24.6 of the labbook 20200707
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/figures", mode: 'copy', pattern: "{*.png}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    // publishDir "${out_path}/reports", mode: 'copy', pattern: "{plot_coverage_report.txt}", overwrite: false // 
    cache 'true'

    input:
    val fastq_name
    path cov // warning: 3 files
    val bef_read_nb
    val after_read_nb
    val ylab
    path cute_file

    output:
    path "plot_${cov.baseName}.png", emit: fig_ch3 // warning: several files
    path "plot_coverage_report.txt", emit: plot_cov_report_ch

    script:
    """
    plot_coverage.R "${cov.baseName}" "${read_nb}" "${xlab}" "${file_name}" "${cute_file}" "plot_coverage_report.txt"
    """

}