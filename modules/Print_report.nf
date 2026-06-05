// The render function only creates the html file with the rmardown
// Inputs:
//      template_rmd: rmardown template file used to create the html (file path is defined in config.nextflow)
//      nb_input: number of fasta sequences in initial input
//      nb_igblast: number of sequences that igblast could analyse
//      nb_unigblast: number of sequences that igblast could not alayse
//      donuts_png: provide links to the donut images needed in the rmd inside the work folder
//      repertoire_png: provide links to the repertoire images needed in the rmd inside the work folder
//      repertoire_constant_ch: names of the constant gene repertoire files to be displayed
//      repertoire_vj_ch: names of the variable gene repertoire files to be displayed
//      itol_subscription: nextflow.config parameter to know if user has paid the subscription to itol automated visualization of trees, process ITOL is only executed if TRUE
//      heavy_chain: to know if the analyzed data is VL or VH, because "Amino acid sequences phylogeny" section in html report is only displayed for VH
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
    path final_warning_ch

    output:
    file "nextflow.config.html"
    file "report.html"
    file "print_report.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    # remove symlink and import folder
    cp ${config_file} config_file.txt
    cp ${template_rmd} report_file.rmd
    if [[ -d "${out_path}/phylo" ]]; then
        cp -r "${out_path}/phylo" .
    else
        mkdir ./phylo
    fi
    if [[ -d "${out_path}/tsv" ]]; then
        cp -r "${out_path}/tsv" .
    else
        mkdir ./tsv
    fi
    if [[ -d "${out_path}/pdf" ]]; then
        cp -r "${out_path}/pdf" .
    else
        mkdir ./pdf
    fi
    cp -r "${out_path}/reports" .
    # end remove symlink and import folder

    Rscript -e '
        warning_collect <- readLines("${final_warning_ch}", warn = FALSE) # one string per line
        rmarkdown::render(
        input = "report_file.rmd",
        output_file = "report.html",
        # list of the variables waiting to be replaced in the rmd file:
        params = list(
            warning_collect = warning_collect
        ),
        # output_dir = ".",
        # intermediates_dir = "./",
        # knit_root_dir = "./",
        run_pandoc = TRUE,
        quiet = TRUE,
        clean = TRUE
        )

        html_here_ok <- TRUE # set to FALSE to rerun the creation of the alignments_viz_html file
        if(html_here_ok == FALSE){
            # dir.create("reports", showWarnings = FALSE, recursive = TRUE)
            rmarkdown::render(input = "alignments_vizu.rmd", output_file = "alignments_viz.html", run_pandoc = TRUE, quiet = TRUE, clean = TRUE)
        }
    ' |& tee -a print_report.log
    """
}
