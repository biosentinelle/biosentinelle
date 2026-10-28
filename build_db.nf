nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     build_db.nf of biosentinelle                                    ##
##                                                                     ##
##     Builds the strain identity database. Run it when the set of     ##
##     reference strains changes, not for every sample.                ##
##                                                                     ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################
*/

include {AssemblyVariants; ReadsIndex; ReadsAlign; ReadsCall; StrainConsensus; PopulationPanel} from './modules/BuildPanel.nf'

/*
Input: the strain catalogue written by bin/build_strain_catalog.py, one row per
strain, with the source its variants come from:

    source_type   source
    matrix        sample name in population_matrix (1011 Yeast Genomes)
    assembly      path of a whole-genome assembly
    reads         comma separated run accessions, read from reads_dir

The reference genome is the coordinate system of the whole database: changing
it invalidates every panel.
*/

workflow {

    print("\n\nINITIATION TIME: ${workflow.start}")
    print("\n\nDATABASE DIRECTORY: ${db_path}/${organism_dir}")
    print("REFERENCE: ${panel_reference}")
    print("K-MER SIZE: ${panel_kmer_size}\n\n")

    reference = file("${panel_reference}", checkIfExists: true)
    matrix = file("${population_matrix}", checkIfExists: true)
    catalogue = file("${strain_table}", checkIfExists: true)

    rows_ch = Channel.fromPath(catalogue).splitCsv(header: true, sep: '\t')

    // assemblies
    assemblies_ch = rows_ch
        .filter{ row -> row.source_type == 'assembly' }
        .map{ row -> tuple(row.strain, file(row.source, checkIfExists: true)) }
    AssemblyVariants(assemblies_ch, reference)
    assembly_strains_ch = AssemblyVariants.out.assembly_vcf_ch.join(AssemblyVariants.out.callable_ch)

    // reads: one alignment and one call per run, then one consensus per strain
    runs_ch = rows_ch
        .filter{ row -> row.source_type == 'reads' }
        .flatMap{ row -> row.source.tokenize(',').collect{ run ->
            def paired = file("${reads_dir}/${run}_1.fastq.gz")
            def fastq = paired.exists() ?
                [paired, file("${reads_dir}/${run}_2.fastq.gz", checkIfExists: true)] :
                [file("${reads_dir}/${run}.fastq.gz", checkIfExists: true)]
            tuple(row.strain, run, fastq)
        } }
    ploidy_ch = rows_ch
        .filter{ row -> row.source_type == 'reads' }
        .map{ row -> tuple(row.strain, row.ploidy ?: '1') }
    ReadsIndex(reference)
    ReadsAlign(runs_ch, ReadsIndex.out.index_ch, reads_min_depth)
    ReadsCall(ReadsAlign.out.bam_ch, reference)
    per_strain_ch = ReadsCall.out.run_vcf_ch
        .join(ReadsAlign.out.callable_ch, by: [0, 1])
        .groupTuple(by: 0)
        .join(ploidy_ch)
        .map{ strain, runs, vcfs, beds, ploidy -> tuple(strain, ploidy, runs, vcfs, beds) }
    StrainConsensus(per_strain_ch)

    // every single-strain source, then the panel
    singles_ch = assembly_strains_ch
        .mix(StrainConsensus.out.strain_ch)
        .toSortedList{ a, b -> a[0] <=> b[0] }
    single_args_ch = singles_ch.map{ list ->
        list.collect{ "--single ${it[0]}=${it[1].name}=${it[2].name}" }.join(' ')
    }
    single_files_ch = singles_ch.map{ list -> list.collect{ [it[1], it[2]] }.flatten() }

    PopulationPanel(
        catalogue,
        reference,
        matrix,
        single_args_ch,
        single_files_ch,
        panel_kmer_size,
        panel_options
    )

    PopulationPanel.out.panel_sites_ch.subscribe{ path ->
        print("\n\nPANEL WRITTEN: ${db_path}/${organism_dir}/panels/${path.name}\n\n")
    }
}
