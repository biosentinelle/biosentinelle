// Native marker genes: what the strain's own URA3, HIS3, LEU2, TRP1, LYS2,
// ADE2, MET15, CAN1, ARG4 and LEU5 look like (the list is marker_loci.tsv of the
// database). Reads are aligned to the whole reference genome, not to the loci
// alone, so a read from a paralog or a repeat gets a low mapping quality instead
// of landing in a marker. Only the reads over the loci are kept afterwards,
// plus the depth over the mitochondrial genome and the read count of every
// chromosome, for the mtDNA (rho, growth on glycerol) status, and the reads
// that place the tool elements in the genome (ToolPosition).
process MarkerAlign {
    label 'bowtie2' // see the withLabel: bowtie2 in the nextflow config file
    cache 'true'

    input:
    val sample_name
    path fastq
    path reference
    path loci_bed
    path mito_genes // mito_genes.tsv of the database, or NO_MITO
    path tools_fasta // tools.fasta of the database, for the positions of the tool elements

    output:
    tuple path("markers.pileup"), path("markers.depth"), path("markers.sam"), emit: marker_evidence_ch
    tuple path("mito.depth"), path("genome.idxstats"), path("read_length.txt"), emit: mito_evidence_ch
    tuple path("tool_clips.tsv.gz"), path("tool_reads.tsv"), path("tool_mates.tsv"), path("tool_chunks.tsv"), path("read_length.txt"), path("read_layout.txt"), emit: tool_evidence_ch
    path "marker_align_report.txt", emit: marker_align_log_ch

    script:
    // Every fastq file, joined for -U: pairing is not needed here, and the
    // second mate is free depth. Mates keep the same read name, which is how
    // ToolPosition finds the mate of a read of a tool element.
    def files = (fastq instanceof List) ? fastq : [fastq]
    def reads = files.join(',')
    def layout = (files.any{ it.name.startsWith('R1__') } ? 'paired' : 'single') +
        (files.any{ it.name.startsWith('S__') } ? ' pieces' : '') // S__: pieces of long sequences
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

    # Tool elements, for their positions in the genome (locate_tools.py, in
    # ToolPosition). (1) The reads of the elements, aligned end to end on
    # tools.fasta. (2) In the genome alignment, with a mapping quality of 20 or
    # more: the reads soft-clipped by 15 bases or more (chromosome, junction,
    # side, aligned part, clipped bases, name), the clipped bases being looked up in
    # the elements afterwards, and the reads named like a read of an element
    # (its mate). (3) Pieces of 50 bases of the elements, one every 10 bases,
    # aligned to the genome: where the reference itself carries an element.
    bowtie2-build --threads ${task.cpus} ${tools_fasta} tools_ref > /dev/null 2>> marker_align_report.txt
    bowtie2 --sensitive -p ${task.cpus} --no-unal --no-hd -x tools_ref -U ${reads} 2>> marker_align_report.txt \\
        | awk -F '\\t' -v OFS='\\t' '{print \$1, \$2, \$3, \$4, \$6, \$10}' > tool_reads.tsv
    samtools view -F 0x904 -q 20 genome.bam | awk -F '\\t' -v OFS='\\t' -v names=tool_reads.tsv -v min_clip=15 '
        BEGIN {
            while ((getline line < names) > 0) { split(line, f, "\\t"); n = f[1]; sub(/\\/[12]\$/, "", n); keep[n] = 1 }
            close(names)
        }
        {
            n = \$1; sub(/\\/[12]\$/, "", n)
            if (n in keep) print \$1, \$2, \$3, \$4, \$6, \$10 > "tool_mates.tsv"
            if (\$6 !~ /S/) next
            lead = 0; trail = 0
            if (match(\$6, /^[0-9]+S/)) lead = substr(\$6, 1, RLENGTH - 1) + 0
            if (match(\$6, /[0-9]+S\$/)) trail = substr(\$6, RSTART, RLENGTH - 1) + 0
            if (lead < min_clip && trail < min_clip) next
            span = 0; cigar = \$6
            while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
                op = substr(cigar, RLENGTH, 1); len = substr(cigar, 1, RLENGTH - 1) + 0
                if (op ~ /[MDN=X]/) span += len
                cigar = substr(cigar, RLENGTH + 1)
            }
            start = \$4 - 1; end = start + span
            if (lead >= min_clip) print \$3, start, "R", start, end, substr(\$10, 1, lead), \$1 > "tool_clips.tsv"
            if (trail >= min_clip) print \$3, end, "L", start, end, substr(\$10, length(\$10) - trail + 1), \$1 > "tool_clips.tsv"
        }'
    touch tool_mates.tsv tool_clips.tsv
    gzip -f tool_clips.tsv
    awk 'function cut(   i, last) {
            last = 0
            for (i = 1; i + 49 <= length(seq); i += 10) { print ">" name ":" i; print substr(seq, i, 50); last = i }
            if (last > 0 && last + 49 < length(seq)) { print ">" name ":" (length(seq) - 49); print substr(seq, length(seq) - 49, 50) }
        }
        /^>/ { if (seq != "") cut(); name = substr(\$1, 2); seq = ""; next }
        { seq = seq toupper(\$0) }
        END { if (seq != "") cut() }' ${tools_fasta} > tool_chunks.fa
    if [[ -s tool_chunks.fa ]]; then
        bowtie2 --local -f -p ${task.cpus} --no-unal --no-hd -x genome_ref -U tool_chunks.fa 2>> marker_align_report.txt \\
            | awk -F '\\t' -v OFS='\\t' '{ s = "NA"; for (i = 12; i <= NF; i++) if (\$i ~ /^AS:i:/) s = substr(\$i, 6); print \$1, \$2, \$3, \$4, \$5, \$6, s }' > tool_chunks.tsv
    else
        touch tool_chunks.tsv
    fi
    echo "${layout}" > read_layout.txt
    echo "reads aligned to the tool elements: \$(wc -l < tool_reads.tsv), soft-clipped genome reads: \$(zcat tool_clips.tsv.gz | wc -l), element pieces on the genome: \$(wc -l < tool_chunks.tsv)" >> marker_align_report.txt
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
