nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     main.nf of biosentinelle                                        ##
##                                                                     ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################
*/

//////// Processes


//////// End Processes


//////// Checking file


//////// end Checking file


//////// Modules

include {CheckVariables} from './conf/CheckVariables.nf'
include {reportEmptyProcess; copyLogFile} from './modules/Functions.nf'
include {Unzip} from './modules/Unzip.nf'
include {WorkflowParam} from './modules/WorkflowParam.nf'
include {Fastq_dump} from './modules/Fastq_dump.nf'
include {Bowtie2} from './modules/Bowtie2.nf'
include {Print_warnings} from './modules/Print_warnings.nf'
include {Print_report} from './modules/Print_report.nf'
include {Backup} from './modules/Backup.nf'

//////// end Modules


//////// Workflow


workflow {

    print("\n\nINITIATION TIME: ${workflow.start}")

    //////// Options of nextflow run

    // --modules (it is just for the process WorkflowParam)
    params.modules = "" // if --module is used, this default value will be overridden
    // end --modules (it is just for the process WorkflowParam)

    //////// end Options of nextflow run


    //////// Variables

    modules = params.modules // remove the dot -> can be used in bash scripts
    config_file = workflow.configFiles[0] // better to use this than config_file = file("${projectDir}/nextflow.config") because the latter is not good if -c option of nextflow run is used
    log_file = file("${launchDir}/.nextflow.log")

    //////// end Variables


    //////// inititation

    print("\n\nRESULT DIRECTORY: ${out_path}")
    if("${system_exec}" == "slurm"){
        print("    queue: ${slurm_queue}")
        print("    qos: ${slurm_qos}")
    }
    if("${system_exec}" != "local"){
        print("    add_options: ${add_options}")
    }
    warning_ch = Channel.empty() // to collect all the warnings
    warn = "\n\nWARNING:\nCURRENTLY FOR Saccharomyces cerevisiae."
    print(warn)
    warning_ch = warning_ch.mix(Channel.value(warn))
    print("\n\n")

    //////// end inititation


    //////// Variable modification

    fastq_name = file("${fastq_path}").baseName
    ref_name = file("${ref_path}").baseName

    //////// end Variable modification


    //////// Channels

    // fs_ch define below because can be a .zip file
    // warning_ch = Channel.empty() // already set above
    // for print_report
    // end for print_report


    //////// end Channels


    //////// files import

    template_rmd = file(template_rmd_path)

    //////// end files import


    //////// Main

    CheckVariables()

    if(fastq_path =~ /.*\.zip$/){
        Unzip( // warning: it is a process defined above
            Channel.fromPath(fastq_path),
            fastq_path
        ) 
        fastq = Unzip.out.unzip_ch.flatten()
    }else if(fastq_path =~ /SRR.*/){
        Fastq_dump(
            fastq_path
        )
        fastq = Fastq_dump.out.fastq_dump_ch.flatten()
        copyLogFile('fastq_dump_report.log', Fastq_dump.out.fastq_dump_log_ch, out_path)
    }else{
        fastq = Channel.fromPath("${fastq_path}", checkIfExists: false) // in channel because many files 
    }

    if(ref_path =~ /.*\.zip$/){
        Unzip( // warning: it is a process defined above
            Channel.fromPath(ref_path),
            ref_path
        ) 
        ref = Unzip.out.unzip_ch.flatten()
    }else{
        ref = Channel.fromPath("${ref_path}", checkIfExists: false) // in channel because many files 
    }


    Bowtie2(
        fastq_name,
        ref_name,
        fastq,
        ref
    )
    copyLogFile('bowtie2_report.log', Bowtie2.out.bowtie2_log_ch, out_path)

    Print_warnings(
        warning_ch.ifEmpty{''}.collectFile(name: "warnings_collect.txt")
    )

/*
    Print_report(
        config_file, // from parameter
        template_rmd, // from parameter
        Print_warnings.out.final_warning_ch // just so that print_report wait for all warnings // warning_ch.collect().map{it.join('\n\n')}.ifEmpty{'EMPTY'} // concatenate all warnings into a single string // finally, the gathered string is very loong. I prefer to use a file added in /reports/
    )


    Backup(
        config_file, 
        log_file
    )
*/

}

    //////// end Main


//////// end Processes