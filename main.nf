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
include {Split_fasta} from './modules/Split_fasta.nf'
include {Fastq_dump} from './modules/Fastq_dump.nf'
include {Bowtie2} from './modules/Bowtie2.nf'
include {VariantCall} from './modules/VariantCall.nf'
include {CountSnps} from './modules/CountSnps.nf'
include {Print_snp_count} from './modules/Print_snp_count.nf'
include {Print_warnings} from './modules/Print_warnings.nf'
include {Print_report} from './modules/Print_report.nf'
include {Q20} from './modules/Q20.nf'
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

    // Note: fastq_name will be set after file extraction (for SRA) or from the fastq_path directly (for local/zip files)
    // This is handled below in the channels section
    ref_name = file("${ref_path}").baseName

    //////// end Variable modification


    //////// Channels

    // ref_ch define below because can be a .zip file
    // warning_ch = Channel.empty() // already set above
    // for print_report
    // end for print_report


    //////// end Channels


    //////// files import

    template_rmd = file(template_rmd_path)

    //////// end files import


    //////// Main

    CheckVariables()

    WorkflowParam(
        modules
    )

    file("${out_path}/tsv").mkdirs()


    if(fastq_path =~ /.*\.zip$/){
        Unzip( // warning: it is a process defined above
            Channel.fromPath(fastq_path),
            fastq_path
        ) 
        fastq = Unzip.out.unzip_ch.flatten()
        fastq_name = file("${fastq_path}").baseName
    }else if(fastq_path =~ /^SRR.*/){
        Fastq_dump(
            fastq_path
        )
        fastq = Fastq_dump.out.fastq_dump_ch.flatten()
        fastq_name = file("${fastq_path}").baseName  // Use SRA accession as the base name
        copyLogFile('fastq_dump_report.log', Fastq_dump.out.fastq_dump_log_ch, out_path)
    }else{
        fastq = Channel.fromPath("${fastq_path}", checkIfExists: false) // in channel because many files
        fastq_name = file("${fastq_path}").baseName
    }
    // is the path a dir or a single file ?
    fastq.branch {
            dir: it.isDirectory()
            file: true
        }.set { branched }
    // Handle directories: list contents
    fastq_ch_from_dir = branched.dir.flatMap { it.listFiles() } // is it is a dir, then recover all the files
    // Handle files: pass through
    fastq_ch_from_file = branched.file
    // Merge back
    fastq_ch = fastq_ch_from_dir.mix(fastq_ch_from_file)
    // end is the path a dir or a single file ?
    fastq_ch.toList().branch {
            single: it.size() == 1
                return it[0]
            multiple: true
                return it
        }.set { branched }




    if(ref_path =~ /.*\.zip$/){
        Unzip( // warning: it is a process defined above
            Channel.fromPath(ref_path),
            ref_path
        ) 
        dir_ch = Unzip.out.unzip_ch.flatten()
    }else{
        dir_ch = Channel.fromPath("${ref_path}", checkIfExists: false) // in channel because many files 
    }
    // is the path a dir or a single file ?
    dir_ch.branch {
            dir: it.isDirectory()
            file: true
        }.set { branched }
    // Handle directories: list contents
    ref_ch_from_dir = branched.dir.flatMap { it.listFiles() } // is it is a dir, then recover all the files
    // Handle files: pass through
    ref_ch_from_file = branched.file
    // Merge back
    ref_ch = ref_ch_from_dir.mix(ref_ch_from_file)
    // end is the path a dir or a single file ?
    ref_ch.toList().branch {
            single: it.size() == 1
                return it[0]
            multiple: true
                return it
        }.set { branched }


    Split_fasta(branched.single)
    ref_ch2 = Split_fasta.out.split_fasta_ch.mix(branched.multiple.flatten()).flatten()

    nb_input = ref_ch2.count()
    ref_ch3 = ref_ch2.map{file -> name = file.baseName ; tuple(file, name)}

   // Align fastq against each reference (parallelization)
    Bowtie2(
        fastq_name,
        fastq_ch.first(),
        ref_ch3
    )
    copyLogFile('bowtie2_report.log', Bowtie2.out.bowtie2_log_ch, out_path)

    // Create bowtie2 input by combining reference info with fastq files
    // For each reference, process against all fastq files
    
    // Create one job per reference that processes all fastq files

    Q20(
        Bowtie2.out.bowtie2_ch
    )
    copyLogFile('q20_report.log', Q20.out.q20_report_ch, out_path)

 
    // Call variants for each reference
    VariantCall(
        Bowtie2.out.bowtie2_ch
    )
    
    // Count SNPs for each reference
    CountSnps(
        VariantCall.out.vcf_ch
    )
    // Collect all SNP counts and find best reference
    Print_snp_count(
        CountSnps.out.snp_count_ch
        .toSortedList { a, b -> a[1] <=> b[1] }   // ascending
        .map { list ->
            def lines = ["Name\tNb_of_variants"] + list.collect { "${it[0]}\t${it[1]}" }
            lines.join("\n") + "\n"
        }
        .collectFile(name: 'snp_counts.tsv')
    )

    Print_warnings(
        warning_ch.ifEmpty{''}.collectFile(name: "warnings_collect.txt")
    )


    Print_report(
        config_file, // from parameter
        template_rmd, // from parameter
        nb_input, // mandatory
        Print_snp_count.out.final_tsv_ch, // just so that print_report wait for all tsv
        Print_warnings.out.final_warning_ch // just so that print_report wait for all warnings // warning_ch.collect().map{it.join('\n\n')}.ifEmpty{'EMPTY'} // concatenate all warnings into a single string // finally, the gathered string is very loong. I prefer to use a file added in /reports/
    )


    Backup(
        config_file, 
        log_file
    )


}

    //////// end Main


//////// end Processes