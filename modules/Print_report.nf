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
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{alignments_viz.html}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{print_report.log}", overwrite: false
    cache 'true'

    input:
    path config_file // to have the file in the work dir
    path template_rmd // to have the file in the work dir
    path alignments_viz_rmd // to have the file in the work dir
    path alignments_viz_html // to have the file in the work dir
    val igblast_organism
    val igblast_loci
    val igblast_B_heavy_chain
    val igblast_B_lambda_chain
    val igblast_B_kappa_chain
    val igblast_T_alpha_chain
    val igblast_T_beta_chain
    val igblast_T_gamma_chain
    val igblast_T_delta_chain
    val clone_strategy
    val clone_model
    val clone_normalize
    val clone_distance
    val clone_germline_kind
    val clone_mut_obs_seq
    val clone_mut_germ_seq
    val clone_mut_regionDefinition
    val align_clone_nb
    val align_soft
    val align_seq
    val align_abalign_options
    val align_mafft_all_options
    val align_mafft_clonal_options
    val nb_input
    val nb_igblast
    val nb_unigblast
    val nb_productive
    val nb_unproductive
    val nb_wanted
    val nb_unwanted
    val nb_dist_ignored
    val nb_clone_tot
    val nb_unclone_tot
    val nb_clone_unassignment
    val nb_clone_ungermline
    path distance_hist_ch // to have the folder in the work dir
    path donuts_png // to have the folder in the work dir
    path repertoire_png // to have the folder in the work dir
    val repertoire_constant_ch
    val repertoire_vj_ch
    val itol_subscription
    path final_warning_ch

    output:
    file "nextflow.config.html"
    file "report.html"
    file "alignments_viz.html"
    file "print_report.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    # remove symlink and import folder
    cp ${config_file} config_file.txt
    cp ${template_rmd} report_file.rmd
    cp ${alignments_viz_rmd} alignments_vizu.rmd
    cp ${alignments_viz_html} alignments_viz.html
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
    cp -r "${projectDir}/bin/doc_images" .
    cp -r "${out_path}/reports" .
    # end remove symlink and import folder

    Rscript -e '
        if(as.integer("${nb_wanted}") > 0){
            # Find the constant and vj repertoires to be displayed in the html file (file names differ depending on light/heavy chain)
            constant_files <- as.character("${repertoire_constant_ch}")
            vj_files <- as.character("${repertoire_vj_ch}")
            cleaned_repertoire_constant <- gsub("^\\\\[\\\\[|\\\\]\\\\]\$", "", constant_files)
            cleaned_repertoire_vj <- gsub("^\\\\[\\\\[|\\\\]\\\\]\$", "", vj_files)
            constant_paths <- strsplit(cleaned_repertoire_constant, ",\\\\s*")[[1]]
            vj_paths <- strsplit(cleaned_repertoire_vj, ",\\\\s*")[[1]]
            constant_names <- basename(constant_paths)
            vj_names <- basename(vj_paths)
            # Verification that the resulting file names are as expected
            constant_rep <- constant_names[grepl("^(IG|TR).C_.*gene_non-zero\\\\.png\$", constant_names)]
            vj_rep <- vj_names[grepl("^rep_gene_(IG|TR).V_.*non-zero\\\\.png\$", vj_names)]
            if(length(constant_rep) == 0 || length(vj_rep) == 0){
                stop(paste0("\\n\\n========\\n\\nERROR IN print_report PROCESS\\n\\nTHE REPERTOIRE PNG FILES TO BE DISPLAYED WERE NOT FOUND\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n"), call. = FALSE)
            }
        }else{
            constant_rep <- ""
            vj_rep <- ""
        }
        # end Find the constant and vj repertoires to be displayed in the html file (file names differ depending on light/heavy chain)
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
            igblast_organism = "${igblast_organism}",
            igblast_loci = "${igblast_loci}",
            igblast_B_heavy_chain = ${igblast_B_heavy_chain},
            igblast_B_lambda_chain = ${igblast_B_lambda_chain},
            igblast_B_kappa_chain = ${igblast_B_kappa_chain},
            igblast_T_alpha_chain = ${igblast_T_alpha_chain},
            igblast_T_beta_chain = ${igblast_T_beta_chain},
            igblast_T_gamma_chain = ${igblast_T_gamma_chain},
            igblast_T_delta_chain = ${igblast_T_delta_chain},
            clone_strategy = "${clone_strategy}",
            clone_model = "${clone_model}",
            clone_normalize = "${clone_normalize}",
            clone_distance = "${clone_distance}",
            clone_germline_kind = "${clone_germline_kind}",
            clone_mut_obs_seq = "${clone_mut_obs_seq}",
            clone_mut_germ_seq = "${clone_mut_germ_seq}",
            clone_mut_regionDefinition = "${clone_mut_regionDefinition}",
            align_clone_nb = "${align_clone_nb}",
            align_soft = "${align_soft}",
            align_seq = "${align_seq}",
            align_abalign_options = "${align_abalign_options}",
            align_mafft_all_options = "${align_mafft_all_options}",
            align_mafft_clonal_options = "${align_mafft_clonal_options}",
            nb_input = ${nb_input},
            nb_igblast = ${nb_igblast}, 
            nb_unigblast = ${nb_unigblast},
            nb_productive = nb_productive, 
            nb_unproductive = nb_unproductive, 
            nb_wanted =  nb_wanted, 
            nb_unwanted =  nb_unwanted, 
            nb_dist_ignored = nb_dist_ignored, 
            nb_clone_tot = nb_clone_tot, 
            nb_unclone_tot = nb_unclone_tot, 
            nb_clone_unassignment = nb_clone_unassignment, 
            nb_clone_ungermline = nb_clone_ungermline, 
            constant_rep = constant_rep,
            vj_rep = vj_rep,
            align_soft = "${align_soft}",
            itol_subscription = ${itol_subscription},
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
