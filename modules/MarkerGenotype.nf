// Native marker genes: what the strain's own URA3, HIS3, LEU2, TRP1, LYS2,
// ADE2, MET15, CAN1, ARG4 and LEU5 look like (the list is marker_loci.tsv of the
// database). Reads are aligned to the whole reference genome, not to the loci
// alone, so a read from a paralog or a repeat gets a low mapping quality instead
// of landing in a marker. Only the reads over the loci are kept afterwards,
// plus the depth over the mitochondrial genome and the read count of every
// chromosome, for the mtDNA (rho, growth on glycerol) status.
process MarkerAlign {
    label 'bowtie2' // see the withLabel: bowtie2 in the nextflow config file
    cache 'true'

    input:
    val sample_name
    path fastq
    path reference
    path loci_bed
    path mito_genes // mito_genes.tsv of the database, or NO_MITO

    output:
    tuple path("markers.pileup"), path("markers.depth"), path("markers.sam"), emit: marker_evidence_ch
    tuple path("mito.depth"), path("genome.idxstats"), path("read_length.txt"), emit: mito_evidence_ch
    path "marker_align_report.txt", emit: marker_align_log_ch

    script:
    // Every fastq file, joined for -U: pairing is not needed here, and the
    // second mate is free depth.
    def reads = (fastq instanceof List) ? fastq.join(',') : fastq
    """
    echo -e "\\n\\n<br /><br />\\n\\n### NATIVE MARKER ALIGNMENT for ${sample_name}\\n\\n" > marker_align_report.txt
    bowtie2-build --threads ${task.cpus} ${reference} genome_ref > /dev/null 2>> marker_align_report.txt

    # Local mode on purpose: reads crossing a junction with foreign sequence (a
    # Ty1, a cassette) are soft-clipped at the junction instead of being dropped,
    # and those clips are the evidence for an insertion.
    bowtie2 --sensitive-local -p ${task.cpus} -x genome_ref -U ${reads} 2>> marker_align_report.txt \\
        | samtools sort -@ ${task.cpus} -o genome.bam -
    samtools index genome.bam
    samtools view -b -L ${loci_bed} genome.bam > markers.bam
    samtools index markers.bam

    # mtDNA: reads per chromosome (nuclear depth), the depth of every base of the
    # mitochondrial chromosome, all mapping qualities, and the mean read length.
    samtools idxstats genome.bam > genome.idxstats
    samtools stats genome.bam | awk -F '\\t' '\$1 == "SN" && \$2 == "average length:" {print \$3}' > read_length.txt
    if [[ "${mito_genes}" != "NO_MITO" ]]; then
        mito_chrom=\$(awk -F '\\t' 'NR == 2 {print \$4}' ${mito_genes})
        samtools depth -a -Q 0 -r "\${mito_chrom}" genome.bam > mito.depth
    else
        touch mito.depth
    fi
    rm genome.bam genome.bam.bai

    # Calls use well mapped reads only. Depth counts every read, so that a region
    # that is merely repetitive is not mistaken for a deletion.
    samtools mpileup -B -q 20 -Q 13 -d 10000 -l ${loci_bed} -f ${reference} markers.bam > markers.pileup 2>> marker_align_report.txt
    samtools depth -a -b ${loci_bed} -Q 0 markers.bam > markers.depth
    samtools view -q 1 markers.bam > markers.sam
    echo "reads over the marker loci: \$(samtools view -c markers.bam)" >> marker_align_report.txt
    """
}

// Call each marker's allele and judge it by the protein it encodes, against the
// allele catalogue built from the strain assemblies.
process MarkerGenotype {
    label 'python' // see the withLabel: python in the nextflow config file
    publishDir "${out_path}/tsv", mode: 'copy', pattern: "marker_genotype.{json,tsv}", overwrite: false
    cache 'true'

    input:
    val sample_name
    tuple path(pileup), path(depth), path(sam)
    path loci
    path loci_fasta
    path alleles
    path documented
    path repeats
    path tools_fasta
    path tools_manifest
    val min_depth
    tuple path(mito_depth), path(idxstats), path(read_length)
    path mito_genes // mito_genes.tsv of the database, or NO_MITO
    val seq_kind
    val seq_molecule

    output:
    path "marker_genotype.json", emit: marker_json_ch
    path "marker_genotype.tsv", emit: marker_tsv_ch
    path "marker_genotype_report.txt", emit: marker_genotype_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### NATIVE MARKER GENOTYPE for ${sample_name}\\n\\n" > marker_genotype_report.txt
    # mtDNA status (rho+, rho-, rho0), hence growth on glycerol.
    mito_arg=""
    if [[ "${mito_genes}" != "NO_MITO" ]]; then
        mito_status.py \\
            --genes ${mito_genes} \\
            --depth ${mito_depth} \\
            --idxstats ${idxstats} \\
            --read-length \$(cat ${read_length}) \\
            --seq-kind ${seq_kind} \\
            --seq-molecule ${seq_molecule} \\
            --min-depth ${min_depth} \\
            --out mito_status.json |& tee -a marker_genotype_report.txt
        mito_arg="--mito mito_status.json"
    fi
    genotype_markers.py \\
        --loci ${loci} \\
        --loci-fasta ${loci_fasta} \\
        --alleles ${alleles} \\
        --documented ${documented} \\
        --repeats ${repeats} \\
        --pileup ${pileup} \\
        --depth ${depth} \\
        --sam ${sam} \\
        --tools ${tools_fasta} \\
        --tools-manifest ${tools_manifest} \\
        --min-depth ${min_depth} \\
        --out marker_genotype.json \\
        --tsv marker_genotype.tsv \${mito_arg} |& tee -a marker_genotype_report.txt
    """
}
