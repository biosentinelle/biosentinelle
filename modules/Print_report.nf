// The render function only creates the html file with the rmardown
// Inputs:
//      template_rmd: rmardown template file used to create the html (file path is defined in config.nextflow)
// Outputs:
//      "report.html": finalized html report for a specific run
//      "print_report.log": will contain any error or warning messages produced by rmardown::render

process Print_report{
    label 'r_ig_clustering'

    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{nextflow.config.html}", overwrite: false
    publishDir path: "${out_path}", mode: 'copy', pattern: "{report.html}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{print_report.log}", overwrite: false
    cache 'true'

    input:
    path config_file // to have the file in the work dir
    path template_rmd // to have the file in the work dir
    val nb_input
    path final_tsv_ch
    path final_tsv_ch
    path final_warning_ch

    output:
    file "nextflow.config.html"
    file "report.html"
    file "print_report.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    cp ${config_file} config_file.txt
    cp ${template_rmd} report_file.rmd
    if [[ -d "${out_path}/tsv" ]]; then
        cp -r "${out_path}/tsv" .
    else
        mkdir ./tsv
    fi
    Rscript -e '
        warning_collect <- readLines("${final_warning_ch}", warn = FALSE) # one string per line
        rmarkdown::render(
            input = "report_file.rmd",
            output_file = "report.html",
            # list of the variables waiting to be replaced in the rmd file:
            params = list(
                nb_input = ${nb_input},
                warning_collect = warning_collect
            ),
            # output_dir = ".",
            # intermediates_dir = "./",
            # knit_root_dir = "./",
            run_pandoc = TRUE,
            quiet = TRUE,
            clean = TRUE
        )

    ' |& tee -a print_report.log
    """
}
