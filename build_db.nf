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

include {AssemblyVariants; BuildSnpPanel; MergePanels} from './modules/BuildPanel.nf'

/*
Input: a two-column TSV, strain name and assembly path.

    strain   assembly
    S288C    /path/to/S288C.fa
    W303     /path/to/W303.fa
    SK1      /path/to/SK1.fa

The reference genome is the coordinate system of the whole database. Every
panel is built against it, so changing it invalidates every panel. The
reference strain itself can be left out of the list: it is the implicit
"carries the reference allele everywhere" member of the panel.
*/

workflow {

    print("\n\nINITIATION TIME: ${workflow.start}")
    print("\n\nDATABASE DIRECTORY: ${db_path}/${organism_dir}")
    print("REFERENCE: ${panel_reference}")
    print("K-MER SIZE: ${panel_kmer_size}\n\n")

    reference = file("${panel_reference}", checkIfExists: true)

    // strain name and assembly file, one row per strain
    assemblies_ch = Channel
        .fromPath("${strain_table}", checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map{ row -> tuple(row.strain, file(row.assembly, checkIfExists: true)) }

    AssemblyVariants(
        assemblies_ch,
        reference
    )

    BuildSnpPanel(
        AssemblyVariants.out.assembly_vcf_ch,
        reference,
        panel_kmer_size,
        panel_max_sites
    )

    // merge_snp_panels.py takes STRAIN=FILE pairs, so the channel is collapsed
    // into the argument string and the matching list of staged files.
    panels_ch = BuildSnpPanel.out.strain_panel_ch.toSortedList{ a, b -> a[0] <=> b[0] }
    callables_ch = AssemblyVariants.out.callable_ch.toSortedList{ a, b -> a[0] <=> b[0] }
    positions_ch = BuildSnpPanel.out.strain_positions_ch.toSortedList{ a, b -> a[0] <=> b[0] }

    panel_args_ch = panels_ch.map{ list ->
        list.collect{ "--panel ${it[0]}=${it[1].name}" }.join(' ')
    }
    callable_args_ch = callables_ch.map{ list ->
        list.collect{ "--callable ${it[0]}=${it[1].name}" }.join(' ')
    }
    position_args_ch = positions_ch.map{ list ->
        list.collect{ "--positions ${it[0]}=${it[1].name}" }.join(' ')
    }

    MergePanels(
        panel_args_ch,
        callable_args_ch,
        position_args_ch,
        reference_strain,
        panel_kmer_size,
        panels_ch.map{ list -> list.collect{ it[1] } },
        callables_ch.map{ list -> list.collect{ it[1] } },
        positions_ch.map{ list -> list.collect{ it[1] } }
    )

    MergePanels.out.panel_sites_ch.subscribe{ path ->
        print("\n\nPANEL WRITTEN: ${db_path}/${organism_dir}/panels/${path.name}\n\n")
    }
}
