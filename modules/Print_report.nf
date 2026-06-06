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
    if [[ -d "${out_path}/tsv" ]]; then
        cp -r "${out_path}/tsv" .
    else
        mkdir ./tsv
    cp -r "${out_path}/reports" .
    # end remove symlink and import folder

    Rscript -e '
        # empty "EMPTY" channel
        nb_productive <- "${nb_productive}"
        nb_unproductive <- "${nb_unproductive}"
        nb_wanted <-  "${nb_wanted}"
        nb_unwanted <-  "${nb_unwanted}"
        nb_dist_ignored <- "${nb_dist_ignored}"
        nb_clone_tot <- "${nb_clone_tot}"
        nb_unclone_tot <- "${nb_unclone_tot}"
        nb_clone_unassignment <- "${nb_clone_unassignment}"
        nb_clone_ungermline <- "${nb_clone_ungermline}"
        warning_collect <- readLines("${final_warning_ch}", warn = FALSE) # one string per line
        if(nb_productive == "EMPTY"){
            nb_productive <- -1
        }else{
            nb_productive <- ${nb_productive}
        }
        if(nb_unproductive == "EMPTY"){
            nb_unproductive <- -1
        }else{
            nb_unproductive <- ${nb_unproductive}
        }
        if(nb_wanted == "EMPTY"){
            nb_wanted <- -1
        }else{
            nb_wanted <- ${nb_wanted}
        }
        if(nb_unwanted == "EMPTY"){
            nb_unwanted <- -1
        }else{
            nb_unwanted <- ${nb_unwanted}
        }
        if(nb_dist_ignored == "EMPTY"){
            nb_dist_ignored <- -1
        }else{
            nb_dist_ignored <- ${nb_dist_ignored}
        }
        if(nb_clone_tot == "EMPTY"){
            nb_clone_tot <- -1
        }else{
            nb_clone_tot <- ${nb_clone_tot}
        }
        if(nb_unclone_tot == "EMPTY"){
            nb_unclone_tot <- -1
        }else{
            nb_unclone_tot <- ${nb_unclone_tot}
        }
        if(nb_clone_unassignment == "EMPTY"){
            nb_clone_unassignment <- -1
        }else{
            nb_clone_unassignment <- ${nb_clone_unassignment}
        }
        if(nb_clone_ungermline == "EMPTY"){
            nb_clone_ungermline <- -1
        }else{
            nb_clone_ungermline <- ${nb_clone_ungermline}
        }
        # empty "EMPTY" channel

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
