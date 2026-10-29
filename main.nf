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
include {VariantCount} from './modules/VariantCount.nf'
include {Sex2tsv} from './modules/Sex2tsv.nf'
include {SexCount} from './modules/SexCount.nf'
include {Print_var_count} from './modules/Print_var_count.nf'
include {Print_sex_count} from './modules/Print_sex_count.nf'
include {Print_warnings} from './modules/Print_warnings.nf'
include {Print_report} from './modules/Print_report.nf'
include {Backup} from './modules/Backup.nf'
include {KmerScan} from './modules/KmerScan.nf'
include {ToolDepth; ToolScore} from './modules/ToolScan.nf'
include {StrainScore; StrainReport} from './modules/StrainID.nf'
include {MarkerAlign; MarkerGenotype} from './modules/MarkerGenotype.nf'

//include {CopyLogFile as CopyLogFile_Closest_germline} from './modules/CopyLogFile.nf'
include {AlignmentToFastq; PrepareReads} from './modules/PrepareReads.nf'
include {ToFasta as ToFasta_ref} from './modules/PrepareReads.nf'
include {ToFasta as ToFasta_feat} from './modules/PrepareReads.nf'
include {ToFasta as ToFasta_sex} from './modules/PrepareReads.nf'
include {Unzip as Unzip_ref} from './modules/Unzip.nf'
include {Unzip as Unzip_feat} from './modules/Unzip.nf'
include {Unzip as Unzip_sex} from './modules/Unzip.nf'
include {Split_fasta as Split_fasta_ref} from './modules/Split_fasta.nf'
include {Split_fasta as Split_fasta_feat} from './modules/Split_fasta.nf'
include {Split_fasta as Split_fasta_sex} from './modules/Split_fasta.nf'
include {Bowtie2 as Bowtie2_ref} from './modules/Bowtie2.nf'
include {Bowtie2 as Bowtie2_feat} from './modules/Bowtie2.nf'
include {Bowtie2 as Bowtie2_sex} from './modules/Bowtie2.nf'
include {Bowtie2 as Bowtie2_tools} from './modules/Bowtie2.nf'
include {Q20 as Q20_ref} from './modules/Q20.nf'
include {Q20 as Q20_feat} from './modules/Q20.nf'
include {Q20 as Q20_sex} from './modules/Q20.nf'
include {VariantCall as VariantCall_ref} from './modules/VariantCall.nf'
include {VariantCall as VariantCall_sex} from './modules/VariantCall.nf'
include {Coverage as Coverage_feat} from './modules/Coverage.nf'
include {Coverage as Coverage_sex} from './modules/Coverage.nf'
include {Plot_coverage as Plot_coverage_feat} from './modules/Plot_coverage.nf'
include {Plot_coverage as Plot_coverage_sex} from './modules/Plot_coverage.nf'
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
    features_name = file("${features_path}").baseName
    sex_name = file("${sex_path}").baseName
    //////// end Variable modification


    //////// Channels

    // ref_ch define below because can be a .zip file
    // warning_ch = Channel.empty() // already set above
    // for print_report
    // end for print_report


    //////// end Channels


    //////// files import

    template_rmd = file(template_rmd_path)
    cute_file = file(cute_path)

    //////// end files import


    //////// Main

    CheckVariables()

    WorkflowParam(
        modules
    )

    file("${out_path}/tsv").mkdirs()

    // sample sequence treatment
    // Any nucleotide format: FASTQ, FASTA, GenBank, EMBL, plain, .gz, .bz2 or
    // .zip, and BAM, SAM, CRAM. Everything ends up as FASTQ through PrepareReads.
    if(fastq_path =~ /^(SRR|ERR|DRR)[0-9]+$/){
        Fastq_dump(
            fastq_path
        )
        fastq = Fastq_dump.out.fastq_dump_ch.flatten()
        copyLogFile('fastq_dump_report.log', Fastq_dump.out.fastq_dump_log_ch, out_path)
    }else{
        fastq = Channel.fromPath("${fastq_path}", checkIfExists: false) // in channel because many files
    }
    fastq_name = file("${fastq_path}").baseName.replaceAll(/\.(fastq|fq|fasta|fa|fna|fas|gb|gbk|gbff|embl|bam|sam|cram)$/, '')
    // is the path a dir or a single file ?
    fastq.branch {
            dir: it.isDirectory()
            file: true
        }.set { branched }
    // Handle directories: list contents, hidden files excluded
    seq_ch = branched.dir.flatMap { it.listFiles().findAll{ f -> ! f.name.startsWith('.') } }.mix(branched.file)
    // end is the path a dir or a single file ?
    seq_ch.branch {
            aln: it.name =~ /(?i)\.(bam|sam|cram)$/
            other: true
        }.set { branched }
    if("${cram_reference}" == "NULL"){
        cram_ref = Channel.of('').collectFile(name: 'NO_CRAM_REF')
    }else{
        cram_ref = Channel.value(file("${cram_reference}", checkIfExists: true))
    }
    AlignmentToFastq(
        branched.aln,
        cram_ref.first()
    )
    copyLogFile('aln_to_fastq_report.log', AlignmentToFastq.out.aln_to_fastq_log_ch, out_path)
    PrepareReads(
        fastq_name,
        branched.other.mix(AlignmentToFastq.out.aln_fastq_ch.flatten()).collect(),
        seq_kind,
        shred_max_length,
        shred_length
    )
    copyLogFile('prepare_reads_report.log', PrepareReads.out.prepare_log_ch, out_path)
    // For the alignments: every file, mates and pieces of long sequences included.
    align_ch = PrepareReads.out.align_ch.collect()
    // For the k-mer scan: every sequence whole, never the pieces.
    scan_ch = PrepareReads.out.align_ch.flatten().filter{ ! it.name.startsWith('S__') }
        .mix(PrepareReads.out.scan_only_ch.flatten())
        .collect()
    // An assembly has a depth of 1 (2 in the overlapping pieces): a read depth
    // threshold would reject every site.
    if("${seq_kind}" == "assembly"){
        kmer_min_depth_used = "1"
        marker_min_depth_used = "1"
    }else{
        kmer_min_depth_used = kmer_min_depth
        marker_min_depth_used = marker_min_depth
    }
    if("${seq_kind}" == "assembly"){
        warn_kind = "\n\nWARNING:\nseq_kind IS \"assembly\": THE COVERAGE PLOTS, READ COUNTS AND VARIANT COUNTS ARE COMPUTED ON PIECES OF THE CONTIGS, NOT ON SEQUENCING READS. THEY DO NOT MEASURE COPY NUMBER.\n"
        print(warn_kind)
        warning_ch = warning_ch.mix(Channel.value(warn_kind))
    }
    if("${seq_molecule}" == "RNA"){
        warn_rna = "\n\nWARNING:\nseq_molecule IS \"RNA\": THE DEPTH OF A GENE FOLLOWS ITS EXPRESSION, NOT ITS COPY NUMBER, IN THE COVERAGE PLOTS AND THE TOOL ELEMENT SCORES. A STRAIN SITE OR A MARKER GENE THAT IS NOT EXPRESSED IS REPORTED AS NOT COVERED.\n"
        print(warn_rna)
        warning_ch = warning_ch.mix(Channel.value(warn_rna))
    }
    // end sample sequence treatment


    // fasta ref treatment
    if(ref_path =~ /.*\.zip$/){
        Unzip_ref( // warning: it is a process defined above
            Channel.fromPath(ref_path),
            ref_path
        ) 
        ref_dir_ch = Unzip_ref.out.unzip_ch.flatten()
    }else{
        ref_dir_ch = Channel.fromPath("${ref_path}", checkIfExists: false) // in channel because many files 
    }
    // is the path a dir or a single file ?
    ref_dir_ch.branch {
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
    // GenBank and EMBL to FASTA
    ToFasta_ref(ref_ch.collect())
    ToFasta_ref.out.fasta_ch.flatten().toList().branch {
            single: it.size() == 1
                return it[0]
            multiple: true
                return it
        }.set { branched }
    Split_fasta_ref(branched.single)
    ref_ch2 = Split_fasta_ref.out.split_fasta_ch.mix(branched.multiple.flatten()).flatten()
    nb_ref = ref_ch2.count()
    ref_ch3 = ref_ch2.map{file -> name = file.baseName ; tuple(file, name)}
    // end fasta ref treatment


    // fasta features treatment
    if(features_path =~ /.*\.zip$/){
        Unzip_feat( // warning: it is a process defined above
            Channel.fromPath(features_path),
            features_path
        ) 
        feat_dir_ch = Unzip_feat.out.unzip_ch.flatten()
    }else{
        feat_dir_ch = Channel.fromPath("${features_path}", checkIfExists: false) // in channel because many files 
    }
    // is the path a dir or a single file ?
    feat_dir_ch.branch {
            dir: it.isDirectory()
            file: true
        }.set { branched }
    // Handle directories: list contents
    feat_ch_from_dir = branched.dir.flatMap { it.listFiles() } // is it is a dir, then recover all the files
    // Handle files: pass through
    feat_ch_from_file = branched.file
    // Merge back
    feat_ch = feat_ch_from_dir.mix(feat_ch_from_file)
    // end is the path a dir or a single file ?
    // GenBank and EMBL to FASTA
    ToFasta_feat(feat_ch.collect())
    ToFasta_feat.out.fasta_ch.flatten().toList().branch {
            single: it.size() == 1
                return it[0]
            multiple: true
                return it
        }.set { branched }
    Split_fasta_feat(branched.single)
    feat_ch2 = Split_fasta_feat.out.split_fasta_ch.mix(branched.multiple.flatten()).flatten()
    nb_feat = feat_ch2.count()
    feat_ch3 = feat_ch2.map{file -> name = file.baseName ; tuple(file, name)}
    // end fasta features treatment


    // fasta sex treatment
    if(sex_path =~ /.*\.zip$/){
        Unzip_sex( // warning: it is a process defined above
            Channel.fromPath(sex_path),
            sex_path
        ) 
        sex_dir_ch = Unzip_sex.out.unzip_ch.flatten()
    }else{
        sex_dir_ch = Channel.fromPath("${sex_path}", checkIfExists: false) // in channel because many files 
    }
    // is the path a dir or a single file ?
    sex_dir_ch.branch {
            dir: it.isDirectory()
            file: true
        }.set { branched }
    // Handle directories: list contents
    sex_ch_from_dir = branched.dir.flatMap { it.listFiles() } // is it is a dir, then recover all the files
    // Handle files: pass through
    sex_ch_from_file = branched.file
    // Merge back
    sex_ch = sex_ch_from_dir.mix(sex_ch_from_file)
    // end is the path a dir or a single file ?
    // GenBank and EMBL to FASTA
    ToFasta_sex(sex_ch.collect())
    ToFasta_sex.out.fasta_ch.flatten().toList().branch {
            single: it.size() == 1
                return it[0]
            multiple: true
                return it
        }.set { branched }
    Split_fasta_sex(branched.single)
    sex_ch2 = Split_fasta_sex.out.split_fasta_ch.mix(branched.multiple.flatten()).flatten()
    nb_sex = sex_ch2.count()
    sex_ch3 = sex_ch2.map{file -> name = file.baseName ; tuple(file, name)}
    // end fasta sex treatment



   // Align fastq against each reference (parallelization)
    Bowtie2_ref(
        fastq_name,
        align_ch,
        ref_ch3
    )
    copyLogFile('bowtie2_ref__report.log', Bowtie2_ref.out.bowtie2_log_ch, out_path)

   // Align fastq against each reference (parallelization)
    Bowtie2_feat(
        fastq_name,
        align_ch,
        feat_ch3
    )
    copyLogFile('bowtie2_feat_report.log', Bowtie2_feat.out.bowtie2_log_ch, out_path)

   // Align fastq against each reference (parallelization)
    Bowtie2_sex(
        fastq_name,
        align_ch,
        sex_ch3
    )
    copyLogFile('bowtie2_sex_report.log', Bowtie2_sex.out.bowtie2_log_ch, out_path)



    Q20_ref(
        Bowtie2_ref.out.bowtie2_ch
    )
    copyLogFile('q20_ref_report.log', Q20_ref.out.q20_report_ch, out_path)

     Q20_feat(
        Bowtie2_feat.out.bowtie2_ch
    )
    copyLogFile('q20_feat_report.log', Q20_feat.out.q20_report_ch, out_path)

     Q20_sex(
        Bowtie2_sex.out.bowtie2_ch
    )
    copyLogFile('q20_sex_report.log', Q20_sex.out.q20_report_ch, out_path)


     Coverage_feat(
        Q20_feat.out.q20_ch
    )
    copyLogFile('coverage_feat_report.log', Coverage_feat.out.cov_report_ch, out_path)


     Plot_coverage_feat(
        Coverage_feat.out.cov_ch,
        Q20_feat.out.q20_bef_read_nb_ch,
        Q20_feat.out.q20_after_read_nb_ch,
        ylab,
        cute_file
    )
    copyLogFile('coverage_feat_report.log', Plot_coverage_feat.out.plot_cov_report_ch, out_path)
    coverage_feat_ch = Plot_coverage_feat.out.cov_tsv_ch.collectFile(name: "coverage_feat.tsv", skip: 1, keepHeader: true)
    coverage_feat_ch.subscribe{it -> it.copyTo("${out_path}/tsv")}


     Coverage_sex(
        Q20_sex.out.q20_ch
    )
    copyLogFile('coverage_sex_report.log', Coverage_sex.out.cov_report_ch, out_path)


     Plot_coverage_sex(
        Coverage_sex.out.cov_ch,
        Q20_sex.out.q20_bef_read_nb_ch,
        Q20_sex.out.q20_after_read_nb_ch,
        ylab,
        cute_file
    )
    copyLogFile('coverage_sex_report.log', Plot_coverage_sex.out.plot_cov_report_ch, out_path)
    coverage_sex_ch = Plot_coverage_sex.out.cov_tsv_ch.collectFile(name: "coverage_sex.tsv", skip: 1, keepHeader: true)
    coverage_sex_ch.subscribe{it -> it.copyTo("${out_path}/tsv")}




    // Call variants for each reference
    VariantCall_ref(
        Q20_ref.out.q20_ch
    )
    copyLogFile('variantCall_ref_report.log', VariantCall_ref.out.variant_call_log_ch, out_path)

    // Call variants for each reference
    VariantCall_sex(
        Q20_sex.out.q20_ch
    )
    copyLogFile('variantCall_sex_report.log', VariantCall_sex.out.variant_call_log_ch, out_path)



    // Count variants for each reference
    VariantCount(
        VariantCall_ref.out.vcf_ch
    )
    Print_var_count(
        VariantCount.out.var_count_ch
        .toSortedList { a, b -> a[1] <=> b[1] }   // ascending
        .map { list ->
            def lines = ["Name\tNb_of_variants"] + list.collect { "${it[0]}\t${it[1]}" }
            lines.join("\n") + "\n"
        }
        .collectFile(name: 'variant_counts.tsv')
    )

/*
    // Count genotypes for each MaT
    Sex2tsv(
        VariantCall_sex.out.vcf_ch
    )
    SexCount(
        Sex2tsv.out.var_count_ch
    )
    Print_sex_count(
        SexCount.out.sex_count_ch
        .map { list ->
            def lines = ["Name\t0:0\t1:0\t1:1"] + list.collect { "${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}" }
            lines.join("\n") + "\n"
        }
        .collectFile(name: 'sex.tsv')
    )
*/


    // Strain identity branch. It is independent of the coverage and variant
    // counting above: the panel carries both allele k-mers, so the sample is
    // never variant called, and the tool elements are scored against the
    // database rather than against whatever was put in features_path.
    if("${strain_id_on}" == "TRUE"){
        db_dir = "${db_path}/${organism_dir}"
        panel_sites = file("${db_dir}/panels/panel_sites.tsv.gz", checkIfExists: false)
        panel_kmers = file("${db_dir}/panels/panel_kmers.tsv.gz", checkIfExists: false)
        panel_strains = file("${db_dir}/panels/panel_strains.tsv", checkIfExists: false)
        mating_kmers = file("${db_dir}/mating/mating_kmers.tsv", checkIfExists: true)
        tools_fasta = file("${db_dir}/tools/tools.fasta", checkIfExists: true)
        tools_manifest = file("${db_dir}/tools/tools.tsv", checkIfExists: true)
        // Reference genome named in the report (e.g. S288C_R64-1-1), from the
        // database's organism.tsv; the organism folder name when it is absent.
        organism_tsv = file("${db_dir}/organism.tsv")
        reference_name = "${organism_dir}"
        if(organism_tsv.exists()){
            ref_row = organism_tsv.readLines().collect{ it.split('\t') }.find{ it[0] == 'reference_name' && it.size() > 1 }
            if(ref_row){ reference_name = ref_row[1] }
        }

        if( ! panel_sites.exists() || ! panel_kmers.exists() || ! panel_strains.exists()){
            warn2 = "\n\nWARNING:\nNO SNP PANEL IN ${db_dir}/panels.\nSTRAIN RANKING IS SKIPPED. RUN build_db.nf FIRST.\n"
            print(warn2)
            warning_ch = warning_ch.mix(Channel.value(warn2))
        }else{
            // Every fastq file, not just the first. The k-mer scan does not use
            // pairing, so the second mate is free depth: on SRR26617475 one mate
            // gives about 11x and both give 23x.
            KmerScan(
                fastq_name,
                scan_ch,
                panel_kmers,
                mating_kmers,
                kmer_max_bases,
                kmer_min_depth_used
            )
            copyLogFile('kmer_scan_report.log', KmerScan.out.kmer_scan_log_ch, out_path)

            Bowtie2_tools(
                fastq_name,
                align_ch,
                Channel.value(tuple(tools_fasta, "tools"))
            )
            copyLogFile('bowtie2_tools_report.log', Bowtie2_tools.out.bowtie2_log_ch, out_path)

            ToolDepth(
                Bowtie2_tools.out.bowtie2_ch
            )
            copyLogFile('tool_depth_report.log', ToolDepth.out.tool_depth_log_ch, out_path)

            ToolScore(
                ToolDepth.out.tool_depth_ch,
                tools_manifest
            )
            copyLogFile('tool_score_report.log', ToolScore.out.tool_score_log_ch, out_path)

            // Native marker genes. The tool scan cannot see them: a wild-type
            // URA3 and a ura3-1 point mutant cover the locus identically.
            marker_loci = file("${db_dir}/markers/marker_loci.tsv", checkIfExists: false)
            marker_alleles = file("${db_dir}/markers/marker_alleles.tsv", checkIfExists: false)
            marker_reference = file("${panel_reference}", checkIfExists: false)
            if( ! marker_loci.exists() || ! marker_alleles.exists() || ! marker_reference.exists()){
                warn3 = "\n\nWARNING:\nNO NATIVE MARKER DATABASE IN ${db_dir}/markers, OR panel_reference NOT FOUND.\nNATIVE MARKER GENOTYPING IS SKIPPED.\n"
                print(warn3)
                warning_ch = warning_ch.mix(Channel.value(warn3))
                marker_json_ch = Channel.of('{}').collectFile(name: 'no_markers.json')
            }else{
                MarkerAlign(
                    fastq_name,
                    align_ch,
                    marker_reference,
                    file("${db_dir}/markers/marker_loci.bed", checkIfExists: true)
                )
                copyLogFile('marker_align_report.log', MarkerAlign.out.marker_align_log_ch, out_path)

                MarkerGenotype(
                    fastq_name,
                    MarkerAlign.out.marker_evidence_ch,
                    marker_loci,
                    file("${db_dir}/markers/marker_loci.fasta", checkIfExists: true),
                    marker_alleles,
                    file("${db_dir}/markers/documented_genotypes.tsv", checkIfExists: true),
                    file("${db_dir}/markers/marker_repeats.tsv", checkIfExists: true),
                    tools_fasta,
                    tools_manifest,
                    marker_min_depth_used
                )
                copyLogFile('marker_genotype_report.log', MarkerGenotype.out.marker_genotype_log_ch, out_path)
                marker_json_ch = MarkerGenotype.out.marker_json_ch
            }

            StrainScore(
                fastq_name,
                KmerScan.out.site_calls_ch,
                KmerScan.out.mating_counts_ch,
                panel_sites,
                panel_strains,
                file("${db_dir}/markers/documented_genotypes.tsv", checkIfExists: true),
                file("${db_dir}/markers/marker_alleles.tsv", checkIfExists: true),
                file("${db_dir}/strains/aliases.tsv", checkIfExists: true),
                tools_manifest,
                ToolScore.out.tool_coverage_ch,
                marker_json_ch,
                organism_name,
                reference_name,
                strain_error_rate,
                min_concordance,
                seq_kind
            )
            copyLogFile('strain_score_report.log', StrainScore.out.strain_score_log_ch, out_path)

            StrainReport(
                fastq_name,
                StrainScore.out.strain_json_ch
            )
        }
    }
    // end Strain identity branch


    Print_warnings(
        warning_ch.ifEmpty{''}.collectFile(name: "warnings_collect.txt")
    )

    // Database files summarized in report.html. Only those that exist are
    // staged; without any (strain identity off, or no database built) a
    // placeholder keeps the input filled and the report says so.
    db_summary_files = []
    if("${strain_id_on}" == "TRUE"){
        db_summary_files = [
            "organism.tsv",
            "panels/panel_sites.tsv.gz",
            "panels/panel_strains.tsv",
            "markers/marker_loci.tsv",
            "markers/marker_loci.fasta",
            "tools/tools.tsv",
            "tools/tools.fasta"
        ].collect{ file("${db_path}/${organism_dir}/${it}") }.findAll{ it.exists() }
    }
    if(db_summary_files.size() == 0){
        db_summary_ch = Channel.of('').collectFile(name: 'no_db.txt')
    }else{
        db_summary_ch = Channel.value(db_summary_files)
    }

    Print_report(
        fastq_name, 
        config_file, // from parameter
        template_rmd, // from parameter
        nb_ref, // mandatory
        nb_feat, // mandatory
        Print_var_count.out.final_tsv_ch, // just so that print_report wait for all tsv
        coverage_feat_ch,
        db_summary_ch, // strain identity database files summarized in the report
        Print_warnings.out.final_warning_ch // just so that print_report wait for all warnings // warning_ch.collect().map{it.join('\n\n')}.ifEmpty{'EMPTY'} // concatenate all warnings into a single string // finally, the gathered string is very loong. I prefer to use a file added in /reports/
    )

    Backup(
        config_file, 
        log_file
    )


}

    //////// end Main


//////// end Processes