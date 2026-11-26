nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     main.nf of lab_strain_watch                                     ##
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
include {WorkflowParam} from './modules/WorkflowParam.nf'
include {Fastq_dump; Fastq_save} from './modules/Fastq_dump.nf'
include {Print_warnings} from './modules/Print_warnings.nf'
include {Print_report} from './modules/Print_report.nf'
include {Backup} from './modules/Backup.nf'
include {KmerScan} from './modules/KmerScan.nf'
include {ToolDepth; ToolScore; ToolPosition} from './modules/ToolScan.nf'
include {StrainScore; StrainReport} from './modules/StrainID.nf'
include {MarkerAlign; MarkerGenotype} from './modules/MarkerGenotype.nf'
include {AlignmentToFastq; PrepareReads} from './modules/PrepareReads.nf'
include {Bowtie2 as Bowtie2_tools} from './modules/Bowtie2.nf'
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
    //////// end Variable modification


    //////// Channels

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

    // sample sequence treatment
    // Any nucleotide format: FASTQ, FASTA, GenBank, EMBL, plain, .gz, .bz2 or
    // .zip, and BAM, SAM, CRAM. Everything ends up as FASTQ through PrepareReads.
    if(fastq_path =~ /^(SRR|ERR|DRR)[0-9]+$/){
        Fastq_dump(
            fastq_path
        )
        fastq = Fastq_dump.out.fastq_dump_ch.flatten()
        copyLogFile('fastq_dump_report.log', Fastq_dump.out.fastq_dump_log_ch, out_path)
        Fastq_save( // gzipped copy in dataset/<accession>/
            fastq_path,
            fastq
        )
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
        warn_kind = "\n\nWARNING:\nseq_kind IS \"assembly\": THE DEPTHS OF THE TOOL ELEMENTS AND OF THE NATIVE MARKER GENES ARE COUNTED ON PIECES OF THE CONTIGS, NOT ON SEQUENCING READS. THEY DO NOT MEASURE COPY NUMBER.\n"
        print(warn_kind)
        warning_ch = warning_ch.mix(Channel.value(warn_kind))
    }
    if("${seq_molecule}" == "RNA"){
        warn_rna = "\n\nWARNING:\nseq_molecule IS \"RNA\": THE DEPTH OF A GENE FOLLOWS ITS EXPRESSION, NOT ITS COPY NUMBER, IN THE TOOL ELEMENT SCORES. A STRAIN SITE OR A MARKER GENE THAT IS NOT EXPRESSED IS REPORTED AS NOT COVERED.\n"
        print(warn_rna)
        warning_ch = warning_ch.mix(Channel.value(warn_rna))
    }
    // end sample sequence treatment


    // Strain identity: the panel carries both allele k-mers, so the sample is
    // never variant called (alignment-free strain and mating type calls), and
    // the reads are aligned only for the tool elements of the database and for
    // the native marker genes, the mitochondrial genome, the 2-micron plasmid
    // and the positions of the tool elements in the genome.
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
        warn2 = "\n\nWARNING:\nNO SNP PANEL IN ${db_dir}/panels.\nSTRAIN RANKING IS SKIPPED. CHECK db_path (A CLONE OF yeast_db) OR RUN ITS build_db.nf FIRST.\n"
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
            warn3 = "\n\nWARNING:\nNO NATIVE MARKER DATABASE IN ${db_dir}/markers, OR panel_reference NOT FOUND.\nNATIVE MARKER GENOTYPING IS SKIPPED, AND THE TOOL ELEMENTS PRESENT ARE NOT PLACED IN THE GENOME.\n"
            print(warn3)
            warning_ch = warning_ch.mix(Channel.value(warn3))
            marker_json_ch = Channel.of('{}').collectFile(name: 'no_markers.json')
            tool_positions_ch = Channel.of('{}').collectFile(name: 'no_tool_positions.json')
        }else{
            // mtDNA status (rho, growth on glycerol): skipped, with a warning,
            // when the database predates mito/mito_genes.tsv.
            mito_genes_file = file("${db_dir}/mito/mito_genes.tsv", checkIfExists: false)
            if(mito_genes_file.exists()){
                mito_genes = Channel.value(mito_genes_file)
            }else{
                warn4 = "\n\nWARNING:\nNO ${db_dir}/mito/mito_genes.tsv.\nTHE MITOCHONDRIAL GENOME (RHO, GROWTH ON GLYCEROL) IS NOT ASSESSED.\n"
                print(warn4)
                warning_ch = warning_ch.mix(Channel.value(warn4))
                mito_genes = Channel.of('').collectFile(name: 'NO_MITO').first()
            }
            // 2-micron plasmid status (cir+, cir0), its sequence being added to
            // the reference genome for the alignment: skipped, with a warning,
            // when the database predates two_micron/.
            two_micron_fasta_file = file("${db_dir}/two_micron/two_micron.fasta", checkIfExists: false)
            two_micron_features_file = file("${db_dir}/two_micron/two_micron_features.tsv", checkIfExists: false)
            if(two_micron_fasta_file.exists() && two_micron_features_file.exists()){
                two_micron_fasta = Channel.value(two_micron_fasta_file)
                two_micron_features = Channel.value(two_micron_features_file)
            }else{
                warn_two_micron = "\n\nWARNING:\nNO ${db_dir}/two_micron/two_micron.fasta OR two_micron_features.tsv.\nTHE 2-MICRON PLASMID (CIR+, CIR0) IS NOT ASSESSED. UPDATE THE CLONE OF yeast_db.\n"
                print(warn_two_micron)
                warning_ch = warning_ch.mix(Channel.value(warn_two_micron))
                two_micron_fasta = Channel.of('').collectFile(name: 'NO_2MICRON').first()
                two_micron_features = two_micron_fasta
            }
            MarkerAlign(
                fastq_name,
                align_ch,
                marker_reference,
                file("${db_dir}/markers/marker_loci.bed", checkIfExists: true),
                mito_genes,
                two_micron_fasta,
                tools_fasta
            )
            copyLogFile('marker_align_report.log', MarkerAlign.out.marker_align_log_ch, out_path)

            // Positions in the genome of the tool elements present, from the
            // same genome alignment: reads crossing a junction between the
            // genome and an element, read pairs with one read in an element,
            // and the native loci of the elements the reference carries.
            ToolPosition(
                MarkerAlign.out.tool_evidence_ch,
                tools_fasta,
                tools_manifest,
                ToolScore.out.tool_coverage_ch,
                marker_reference,
                seq_kind
            )
            copyLogFile('tool_position_report.log', ToolPosition.out.tool_position_log_ch, out_path)
            tool_positions_ch = ToolPosition.out.tool_positions_json_ch

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
                marker_min_depth_used,
                MarkerAlign.out.mito_evidence_ch,
                mito_genes,
                MarkerAlign.out.two_micron_evidence_ch,
                two_micron_features,
                seq_kind,
                seq_molecule
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
            tool_positions_ch,
            organism_name,
            reference_name,
            strain_error_rate,
            min_concordance,
            seq_kind
        )
        copyLogFile('strain_score_report.log', StrainScore.out.strain_score_log_ch, out_path)

        // Tree of the panel strains, opened from the report (button next to
        // "Strains in the panel"): made without it, with a warning, when the
        // database predates panels/panel_tree.nwk.
        panel_tree_file = file("${db_dir}/panels/panel_tree.nwk", checkIfExists: false)
        if(panel_tree_file.exists()){
            panel_tree = Channel.value(panel_tree_file)
        }else{
            warn5 = "\n\nWARNING:\nNO ${db_dir}/panels/panel_tree.nwk.\nTHE STRAIN IDENTITY REPORT HAS NO TREE OF THE PANEL STRAINS. UPDATE THE CLONE OF yeast_db.\n"
            print(warn5)
            warning_ch = warning_ch.mix(Channel.value(warn5))
            panel_tree = Channel.of('').collectFile(name: 'NO_TREE').first()
        }

        StrainReport(
            fastq_name,
            StrainScore.out.strain_json_ch,
            panel_tree,
            panel_strains
        )
    }
    // end Strain identity


    Print_warnings(
        warning_ch.ifEmpty{''}.collectFile(name: "warnings_collect.txt")
    )

    // Database files summarized in report.html. Only those that exist are
    // staged; without any (no database built) a placeholder keeps the input
    // filled and the report says so.
    db_summary_files = [
        "organism.tsv",
        "panels/panel_sites.tsv.gz",
        "panels/panel_strains.tsv",
        "markers/marker_loci.tsv",
        "markers/marker_loci.fasta",
        "tools/tools.tsv",
        "tools/tools.fasta"
    ].collect{ file("${db_path}/${organism_dir}/${it}") }.findAll{ it.exists() }
    if(db_summary_files.size() == 0){
        db_summary_ch = Channel.of('').collectFile(name: 'no_db.txt')
    }else{
        db_summary_ch = Channel.value(db_summary_files)
    }

    Print_report(
        fastq_name, 
        config_file, // from parameter
        template_rmd, // from parameter
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