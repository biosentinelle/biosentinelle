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

    # paftools.js call needs the alignments sorted by reference name and start.
    sort -k6,6 -k8,8n ${strain}.paf > ${strain}.sorted.paf

    # Called twice on purpose. With -f the output is VCF, which is what the
    # panel builder reads. Without -f it is paftools' own format, whose R lines
    # give the reference intervals the assembly actually covered; the VCF does
    # not carry those. Those intervals are what separates "this strain carries
    # the reference allele here" from "this strain has no call here", and
    # merge_snp_panels.py needs the difference to avoid claiming knowledge it
    # does not have.
    paftools.js call -f ${reference} -l 1000 -L 10000 ${strain}.sorted.paf \\
        > ${strain}.vcf 2>> ${strain}_assembly_variants_report.txt

    paftools.js call -l 1000 -L 10000 ${strain}.sorted.paf 2>> ${strain}_assembly_variants_report.txt \\
        | awk -F'\\t' '\$1 == "R" {print \$2"\\t"\$3"\\t"\$4}' \\
        | sort -k1,1 -k2,2n > ${strain}.callable.bed

    if [ ! -s ${strain}.callable.bed ]; then
        echo "ERROR: no callable interval extracted for ${strain}" >> ${strain}_assembly_variants_report.txt
        exit 1
    fi
    echo "callable intervals: \$(wc -l < ${strain}.callable.bed)" >> ${strain}_assembly_variants_report.txt

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

// Index the reference once for every read run.
process ReadsIndex {
    label 'bowtie2' // see the withLabel: bowtie2 in the nextflow config file
    cache 'true'

    input:
    path reference

    output:
    path "reference_index*", emit: index_ch

    script:
    """
    bowtie2-build --threads ${task.cpus} ${reference} reference_index > /dev/null
    """
}

// Align one read run to the reference, and write where it can be read.
process ReadsAlign {
    label 'bowtie2' // see the withLabel: bowtie2 in the nextflow config file
    cache 'true'

    input:
    tuple val(strain), val(run), path(fastq)
    path index
    val min_depth

    output:
    tuple val(strain), val(run), path("${run}.bam"), path("${run}.bam.bai"), emit: bam_ch
    tuple val(strain), val(run), path("${run}.callable.bed"), emit: callable_ch
    path "${run}_align_report.txt", emit: align_log_ch

    script:
    def reads = fastq instanceof List && fastq.size() == 2 ? "-1 ${fastq[0]} -2 ${fastq[1]}" : "-U ${fastq instanceof List ? fastq.join(',') : fastq}"
    """
    echo -e "\n\n### READ ALIGNMENT ${strain} ${run}\n" > ${run}_align_report.txt
    # end-to-end: a soft-clipped read end is where a variant call goes wrong,
    # and these are clean Illumina runs of the strain the reference comes from
    # or a close relative.
    bowtie2 --sensitive -p ${task.cpus} -x reference_index ${reads} 2>> ${run}_align_report.txt \
        | samtools sort -@ 4 -o ${run}.bam -
    samtools index ${run}.bam

    # Callable: positions with at least min_depth reads of mapping quality 20
    # and at most 2.5 times the median depth. A collapsed repeat or a
    # duplication shows as excess depth, and its calls mix the copies.
    samtools depth -a -Q 20 -q 20 ${run}.bam > depth.tsv
    median=\$(cut -f3 depth.tsv | sort -n | awk '{a[NR]=\$1} END{print a[int((NR+1)/2)]}')
    awk -v lo=${min_depth} -v hi=\$(awk -v m=\$median 'BEGIN{print int(2.5*m)}') -F'\t' '
        \$3 >= lo && \$3 <= hi {
            if (\$1 == c && \$2 == e + 1) { e = \$2; next }
            if (c != "") print c"\t"s-1"\t"e
            c = \$1; s = \$2; e = \$2; next
        }
        END { if (c != "") print c"\t"s-1"\t"e }' depth.tsv > ${run}.callable.bed
    rm depth.tsv
    echo "median depth: \$median" >> ${run}_align_report.txt
    echo "callable bases: \$(awk '{s+=\$3-\$2} END{print s}' ${run}.callable.bed)" >> ${run}_align_report.txt
    """
}

// Call one run's variants.
process ReadsCall {
    label 'bcftools' // see the withLabel: bcftools in the nextflow config file
    cache 'true'

    input:
    tuple val(strain), val(run), path(bam), path(bai)
    path reference

    output:
    tuple val(strain), val(run), path("${run}.vcf.gz"), emit: run_vcf_ch

    script:
    """
    bcftools mpileup -f ${reference} -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP -Ou ${bam} \
        | bcftools call -m -v --ploidy 2 -Oz -o ${run}.vcf.gz
    """
}

// One strain from its runs: the variants every run agrees on.
process StrainConsensus {
    label 'python' // see the withLabel: python in the nextflow config file
    publishDir "${db_path}/${organism_dir}/strains/reads", mode: 'copy', overwrite: true
    cache 'true'

    input:
    tuple val(strain), val(ploidy), val(runs), path(vcfs), path(beds)

    output:
    tuple val(strain), path("${strain}.vcf"), path("${strain}.callable.bed"), emit: strain_ch
    path "${strain}_consensus_report.txt", emit: consensus_log_ch

    script:
    def run_args = [runs, vcfs, beds].transpose().collect{ "--run ${it[0]}=${it[1]}=${it[2]}" }.join(' ')
    """
    consensus_strain_calls.py --strain ${strain} --ploidy ${ploidy} ${run_args} \
        --out-vcf ${strain}.vcf --out-bed ${strain}.callable.bed |& tee ${strain}_consensus_report.txt
    """
}

// Every strain at once: the per-strain k-mer panel.
process PopulationPanel {
    label 'python' // see the withLabel: python in the nextflow config file
    publishDir "${db_path}/${organism_dir}/panels", mode: 'copy', overwrite: true
    cache 'true'
    memory '48 GB'

    input:
    path strain_table
    path reference
    path matrix
    val single_args
    path single_files
    val kmer_size
    val panel_options

    output:
    path "panel_sites.tsv.gz", emit: panel_sites_ch
    path "panel_kmers.tsv.gz", emit: panel_kmers_ch
    path "panel_strains.tsv", emit: panel_strains_ch
    path "panel_report.txt", emit: panel_log_ch

    script:
    """
    build_population_panel.py \
        --strains ${strain_table} \
        --reference ${reference} \
        --matrix ${matrix} \
        ${single_args} \
        -k ${kmer_size} \
        ${panel_options} \
        --out-sites panel_sites.tsv \
        --out-kmers panel_kmers.tsv \
        --out-strains panel_strains.tsv |& tee panel_report.txt
    gzip -n panel_sites.tsv panel_kmers.tsv
    """
}
