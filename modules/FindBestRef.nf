// Find best reference (lowest SNP count) and generate summary report
process FindBestRef {
    label 'bash'
    cache 'true'

    input:
    val snp_counts_list  // Will be a list of [ref_name, snp_count] tuples as strings

    output:
    path "best_reference_report.txt", emit: best_ref_report_ch

    script:
    template 'find_best_ref.sh'
}

