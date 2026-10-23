// Database construction. These processes are run by build_db.nf, not by the
// analysis pipeline: they turn a set of strain assemblies into the k-mer SNP
// panel that main.nf then uses read-only.

// Align one strain assembly to the reference and call the fixed differences.
process AssemblyVariants {
    label 'minimap2' // see the withLabel: minimap2 in the nextflow config file
    cache 'true'

    input:
    tuple val(strain), path(assembly)
    path reference

    output:
    tuple val(strain), path("${strain}.vcf"), emit: assembly_vcf_ch
    tuple val(strain), path("${strain}.callable.bed"), emit: callable_ch
    path "${strain}_assembly_variants_report.txt", emit: assembly_variants_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### ASSEMBLY ALIGNMENT ${strain}\\n\\n" > ${strain}_assembly_variants_report.txt

    # asm5 is the right preset for two strains of one species: it expects a few
    # percent divergence at most, which holds even for SK1 against S288C.
    minimap2 -cx asm5 --cs ${reference} ${assembly} > ${strain}.paf 2>> ${strain}_assembly_variants_report.txt

    sort -k6,6 -k8,8n ${strain}.paf > ${strain}.sorted.paf
    paftools.js call -f ${reference} -l 1000 -L 10000 ${strain}.sorted.paf > ${strain}.vcf 2>> ${strain}_assembly_variants_report.txt

    # The regions the assembly actually aligned to. Everything outside is a
    # region where this strain has no call, which is different from a region
    # where it carries the reference allele. merge_snp_panels.py uses this to
    # avoid claiming knowledge it does not have.
    awk -F'\\t' '\$1 ~ /^R/ {print \$2"\\t"\$3"\\t"\$4}' ${strain}.vcf \\
        | sort -k1,1 -k2,2n > ${strain}.callable.bed || touch ${strain}.callable.bed

    echo "variant lines: \$(grep -vc '^#' ${strain}.vcf)" >> ${strain}_assembly_variants_report.txt
    """
}

// Turn one strain's variants into reference and alternative allele k-mers.
process BuildSnpPanel {
    label 'python' // see the withLabel: python in the nextflow config file
    cache 'true'

    input:
    tuple val(strain), path(vcf)
    path reference
    val kmer_size
    val max_sites

    output:
    tuple val(strain), path("${strain}.snp.tsv"), emit: strain_panel_ch
    tuple val(strain), path("${strain}.positions.tsv"), emit: strain_positions_ch
    path "${strain}_panel_report.txt", emit: panel_log_ch

    script:
    """
    build_snp_panel.py \\
        --vcf ${vcf} \\
        --reference ${reference} \\
        --strain ${strain} \\
        --out ${strain}.snp.tsv \\
        --out-positions ${strain}.positions.tsv \\
        -k ${kmer_size} \\
        --max-sites ${max_sites} |& tee ${strain}_panel_report.txt
    """
}

// Merge the per-strain panels into the discriminant site matrix.
process MergePanels {
    label 'python' // see the withLabel: python in the nextflow config file
    publishDir "${db_path}/${organism_dir}/panels", mode: 'copy', overwrite: true
    cache 'true'

    input:
    val panel_args
    val callable_args
    val position_args
    val reference_strain
    val kmer_size
    path panels
    path callables
    path positions

    output:
    path "panel_sites.tsv", emit: panel_sites_ch
    path "panel_kmers.tsv", emit: panel_kmers_ch
    path "merge_panels_report.txt", emit: merge_log_ch

    script:
    """
    merge_snp_panels.py \\
        ${panel_args} \\
        ${callable_args} \\
        ${position_args} \\
        --reference-strain ${reference_strain} \\
        -k ${kmer_size} \\
        --out-sites panel_sites.tsv \\
        --out-kmers panel_kmers.tsv |& tee merge_panels_report.txt
    """
}
