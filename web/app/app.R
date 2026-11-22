#!/usr/bin/Rscript
#########################################################################
##                                                                     ##
##     app.R                                                           ##
##                                                                     ##
##     Web page of lab_strain_watch                                    ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################


################################ Aim

# Simple web form generated from the nextflow.config file of lab_strain_watch:
# one input per parameter of the "Parameters that must be set by the user" section (except those of galaxy/locked_params.txt),
# with the default value and the comment of nextflow.config as help.
# A submission writes the sample files in <LSW_WORK_DIR>/jobs/<key>/sample/, the parameters in <LSW_WORK_DIR>/jobs/<key>/params.json
# and the task in <LSW_WORK_DIR>/todo/<key>.json (with the run accession, if given instead of files).
# The web/daemon/lsw_bioblend.py daemon sends the task to the lab_strain_watch Galaxy tool (galaxy.pasteur.fr -> Maestro cluster).
# The key, given to the user, is the only way to get the results back ("Check my job" tab, or the link in the e-mail).

################################ End Aim


################################ Settings (environment variables)

work_dir <- Sys.getenv("LSW_WORK_DIR", "/srv/lsw_work") # folder shared with the daemon
repo_dir <- Sys.getenv("LSW_REPO_DIR", normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)) # lab_strain_watch clone (nextflow.config and galaxy/)
python <- Sys.getenv("LSW_PYTHON", "python3")
max_upload_mb <- as.numeric(Sys.getenv("LSW_MAX_UPLOAD_MB", "1000"))
retention_days <- Sys.getenv("LSW_RETENTION_DAYS", "15")
options(shiny.maxRequestSize = max_upload_mb * 1024^2)

config_file <- file.path(repo_dir, "nextflow.config")
fill_config <- file.path(repo_dir, "galaxy", "fill_config.py")
locked_file <- file.path(repo_dir, "galaxy", "locked_params.txt")
key_regex <- "^[A-Za-z0-9_-]{16,64}$"
mail_regex <- "^[^@[:space:]]+@[^@[:space:]]+\\.[^@[:space:]]+$"
accession_regex <- "^(SRR|ERR|DRR)[0-9]+$"
sample_regex <- "^[A-Za-z0-9_-]*$"
# sample files read by the pipeline (bin/prepare_reads.py and samtools), and files accepted inside a zip
seq_regex <- "\\.(fastq|fq|fasta|fa|fna|fas|fsa|seq|txt|gb|gbk|gbff|embl|dat)(\\.(gz|bz2))?$"
file_regex <- "(\\.(fastq|fq|fasta|fa|fna|fas|fsa|seq|txt|gb|gbk|gbff|embl|dat)(\\.(gz|bz2))?|\\.(bam|sam|cram|zip))$"
# parameters with a fixed set of values (see their comment in nextflow.config)
choices <- list(seq_kind = c("reads", "assembly"), seq_molecule = c("DNA", "RNA"))

for(i0 in c("todo", "doing", "done", "error", "jobs")){
    dir.create(file.path(work_dir, i0), recursive = TRUE, showWarnings = FALSE)
}

################################ End Settings


################################ Parameters of nextflow.config

# parsed by fill_config.py, the same script that fills the config in Galaxy, so that the page and the tool always agree
spec <- jsonlite::fromJSON(paste(system2(python, c(shQuote(fill_config), "--template", shQuote(config_file), "--spec"), stdout = TRUE), collapse = "\n"), simplifyDataFrame = FALSE)
locked <- trimws(sub("#.*", "", readLines(locked_file)))
locked <- locked[locked != ""]
spec <- spec[ ! sapply(spec, function(x){x$name}) %in% locked]
sections <- unique(sapply(spec, function(x){x$section}))

param_input <- function(p){
    id <- paste0("p_", p$name)
    if(p$default %in% c("TRUE", "FALSE")){
        widget <- shiny::selectInput(id, p$name, choices = c("TRUE", "FALSE"), selected = p$default)
    }else if( ! is.null(choices[[p$name]])){
        widget <- shiny::selectInput(id, p$name, choices = choices[[p$name]], selected = p$default)
    }else{
        widget <- shiny::textInput(id, p$name, value = p$default, width = "100%")
    }
    shiny::div(
        class = "lsw-param",
        widget,
        shiny::tags$details(shiny::tags$summary("help"), shiny::tags$pre(class = "lsw-help", p$help))
    )
}

section_panel <- function(title){
    params <- Filter(function(x){x$section == title}, spec)
    shiny::tags$details(
        open = if(grepl("^Optional", title)){NULL}else{"open"}, # optional sections closed
        class = "lsw-section",
        shiny::tags$summary(shiny::h4(title, style = "display: inline")),
        lapply(params, param_input)
    )
}

################################ End Parameters of nextflow.config


################################ Functions

new_key <- function(){
    # cryptographic random key: it is the only protection of the results
    bytes <- tryCatch({
        con <- file("/dev/urandom", open = "rb", raw = TRUE) # raw = TRUE: device file, not a regular file
        on.exit(close(con))
        readBin(con, "raw", 18)
    }, error = function(e){NULL})
    if(is.null(bytes) || length(bytes) != 18){
        stop("cannot read /dev/urandom to create the job key")
    }
    paste(as.character(bytes), collapse = "")
}

find_task <- function(key){
    # return list(folder, task) or NULL
    if( ! grepl(key_regex, key)){
        return(NULL)
    }
    for(i0 in c("todo", "doing", "done", "error")){
        path <- file.path(work_dir, i0, paste0(key, ".json"))
        if(file.exists(path)){
            task <- tryCatch(jsonlite::fromJSON(path), error = function(e){list(state = "unknown", message = "Status being updated, please retry.")})
            return(list(key = key, folder = i0, task = task))
        }
    }
    NULL
}

safe_names <- function(names){
    # file names kept (the mates are paired by their _1/_2 or _R1/_R2 suffix), made path safe and unique
    names <- gsub("[^A-Za-z0-9._-]", "_", basename(names))
    names <- sub("^[._]+", "", names)
    names[names == ""] <- "sequences"
    dup <- duplicated(names)
    names[dup] <- paste0(seq_along(names)[dup], "_", names[dup])
    names
}

################################ End Functions


################################ UI

ui <- shiny::navbarPage(
    title = "lab_strain_watch",
    id = "tabs",
    header = shiny::tags$style(shiny::HTML("
        .lsw-section {border: 1px solid #ddd; border-radius: 4px; padding: 8px 12px; margin-bottom: 10px;}
        .lsw-param {margin-bottom: 6px;}
        .lsw-help {white-space: pre-wrap; font-size: 12px; max-height: 250px; overflow-y: auto;}
        details > summary {cursor: pointer;}
    ")),
    shiny::tabPanel(
        "Submit",
        value = "submit",
        shiny::fluidRow(shiny::column(
            width = 10, offset = 1,
            shiny::p("Check the identity of a sequenced laboratory yeast strain with the ", shiny::a("lab_strain_watch", href = "https://github.com/lab-strain-watch/lab_strain_watch", target = "_blank"),
                " pipeline, on the Institut Pasteur cluster: most probable strain, mating type, native marker genes, mitochondrial genome and laboratory tool elements.",
                " Fill the form below, which corresponds to the nextflow.config file of the pipeline, and submit.",
                " You will get a key to follow the job and download the results, which are kept ", retention_days, " days."),
            shiny::tags$details(
                open = "open",
                class = "lsw-section",
                shiny::tags$summary(shiny::h4("Sample", style = "display: inline")),
                shiny::fileInput("sample_files", paste0("Sequence file(s) of a single sample (fastq_path parameter, max ", max_upload_mb, " MB in all): FASTQ, FASTA, GenBank, EMBL (plain, .gz or .bz2), BAM, SAM, CRAM, or a zip of FASTQ, FASTA, GenBank or EMBL files. Paired files are recognized by their _1/_2 or _R1/_R2 suffix."), multiple = TRUE, width = "100%"),
                shiny::textInput("accession", "Or a run accession of SRA (SRR...), ENA (ERR...) or DDBJ (DRR...), instead of files"),
                shiny::textInput("sample_name", "Sample name (optional, letters, digits, _ and -; default: name of the first file, or the run accession)"),
                shiny::textInput("mail", "E-mail (optional, to be warned when the job ends)")
            ),
            lapply(sections, section_panel),
            shiny::actionButton("submit", "Submit", class = "btn-primary btn-lg"),
            shiny::br(), shiny::br()
        ))
    ),
    shiny::tabPanel(
        "Check my job",
        value = "check",
        shiny::fluidRow(shiny::column(
            width = 8, offset = 2,
            shiny::textInput("key", "Job key", width = "100%"),
            shiny::actionButton("check", "Check"),
            shiny::hr(),
            shiny::uiOutput("status")
        ))
    )
)

################################ End UI


################################ Server

server <- function(input, output, session){

    # key given in the URL (?key=...), e.g. from the e-mail
    shiny::observe({
        query <- shiny::parseQueryString(session$clientData$url_search)
        if( ! is.null(query$key)){
            shiny::updateTextInput(session, "key", value = query$key)
            shiny::updateNavbarPage(session, "tabs", selected = "check")
        }
    })

    shiny::observeEvent(input$submit, {
        fail <- function(...){
            shiny::showModal(shiny::modalDialog(title = "Submission refused", paste0(...), easyClose = TRUE))
        }
        # sample: files or run accession
        accession <- trimws(input$accession)
        files <- input$sample_files
        if(is.null(files) && accession == ""){
            return(fail("Please give the sequence file(s) of the sample, or a run accession."))
        }
        if( ! is.null(files) && accession != ""){
            return(fail("Please give either sequence files or a run accession, not both."))
        }
        if(accession != "" && ! grepl(accession_regex, accession)){
            return(fail("Invalid run accession: SRR, ERR or DRR followed by digits is expected, e.g., SRR26617475."))
        }
        if( ! is.null(files)){
            bad <- files$name[ ! grepl(file_regex, files$name, ignore.case = TRUE)]
            if(length(bad) > 0){
                return(fail("Not a sequence file (see the accepted formats): ", paste(utils::head(bad, 5), collapse = ", ")))
            }
            for(i1 in which(grepl("\\.zip$", files$name, ignore.case = TRUE))){
                zip_content <- tryCatch(utils::unzip(files$datapath[i1], list = TRUE)$Name, error = function(e){NULL})
                if(length(zip_content) == 0){
                    return(fail(files$name[i1], " is not a valid zip file."))
                }
                not_seq <- zip_content[ ! (grepl(seq_regex, zip_content, ignore.case = TRUE) | grepl("/$", zip_content))]
                if(length(not_seq) > 0){
                    return(fail("The zip files must contain only FASTQ, FASTA, GenBank or EMBL files. Not accepted in ", files$name[i1], ": ", paste(utils::head(not_seq, 5), collapse = ", ")))
                }
            }
        }
        sample_name <- trimws(input$sample_name)
        if( ! grepl(sample_regex, sample_name)){
            return(fail("Invalid sample name: letters, digits, _ and - only."))
        }
        mail <- trimws(input$mail)
        if(mail != "" && ! grepl(mail_regex, mail)){
            return(fail("Invalid e-mail address."))
        }
        # parameters: only those changed from the nextflow.config default
        values <- list()
        for(p in spec){
            value <- input[[paste0("p_", p$name)]]
            if(is.null(value)){ # input not created yet
                next
            }
            value <- trimws(value)
            if( ! identical(value, p$default)){
                values[[p$name]] <- value
            }
        }
        key <- new_key()
        job_dir <- file.path(work_dir, "jobs", key)
        dir.create(job_dir, recursive = TRUE)
        params_file <- file.path(job_dir, "params.json")
        writeLines(jsonlite::toJSON(if(length(values) == 0){structure(list(), names = character(0))}else{values}, auto_unbox = TRUE, pretty = TRUE), params_file)
        # same check as in the Galaxy tool, so that errors are given now and not after the run
        check <- suppressWarnings(system2(python, c(shQuote(fill_config), "--template", shQuote(config_file), "--locked", shQuote(locked_file), "--values", shQuote(params_file), "-o", if(.Platform$OS.type == "windows"){"NUL"}else{"/dev/null"}), stdout = TRUE, stderr = TRUE))
        if( ! is.null(attr(check, "status"))){
            unlink(job_dir, recursive = TRUE)
            return(fail(paste(check, collapse = "\n")))
        }
        if( ! is.null(files)){
            dir.create(file.path(job_dir, "sample"))
            file.copy(files$datapath, file.path(job_dir, "sample", safe_names(files$name)))
        }
        # the task file is written last: the daemon only takes complete jobs
        task <- list(key = key, mail = mail, accession = accession, sample_name = sample_name, submitted = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), state = "queued", message = "The job is waiting to be sent to Galaxy.")
        tmp <- file.path(work_dir, "todo", paste0(key, ".json.tmp"))
        writeLines(jsonlite::toJSON(task, auto_unbox = TRUE, pretty = TRUE), tmp)
        file.rename(tmp, file.path(work_dir, "todo", paste0(key, ".json")))
        shiny::updateTextInput(session, "key", value = key)
        shiny::showModal(shiny::modalDialog(
            title = "Job submitted",
            shiny::p("Keep this key to follow your job and download the results:"),
            shiny::tags$pre(key),
            shiny::p("You can close this page and come back later with the key (Check my job tab)."),
            easyClose = TRUE
        ))
    })

    current <- shiny::eventReactive(input$check, {
        find_task(trimws(input$key))
    })

    output$status <- shiny::renderUI({
        res <- current()
        if(is.null(res)){
            return(shiny::p("Unknown key."))
        }
        items <- list(
            shiny::h4("State: ", res$task$state),
            shiny::p(res$task$message),
            shiny::p("Submitted: ", res$task$submitted)
        )
        if(res$folder == "done"){
            items <- c(items, list(
                shiny::downloadButton("dl_strain", "Strain identity report (.html)"), " ",
                shiny::downloadButton("dl_report", "Coverage report (.html)"), " ",
                shiny::downloadButton("dl_results", "All the results (.zip)")
            ))
        }else if(res$folder == "error" && file.exists(file.path(work_dir, "jobs", res$key, "error.txt"))){
            items <- c(items, list(shiny::downloadButton("dl_error", "Download the error log")))
        }else if(res$folder != "error"){
            items <- c(items, list(shiny::p(shiny::em("Click Check again to refresh."))))
        }
        do.call(shiny::tagList, items)
    })

    job_file <- function(name){
        shiny::downloadHandler(
            filename = function(){paste0("lab_strain_watch_", name)},
            content = function(file){
                res <- current()
                shiny::req(res)
                path <- file.path(work_dir, "jobs", res$key, name)
                shiny::req(file.exists(path)) # e.g., no strain identity report when the database has no strain panel
                file.copy(path, file)
            }
        )
    }
    output$dl_strain <- job_file("strain_identity_report.html")
    output$dl_report <- job_file("report.html")
    output$dl_results <- job_file("results.zip")
    output$dl_error <- job_file("error.txt")
}

################################ End Server

shiny::shinyApp(ui = ui, server = server)
