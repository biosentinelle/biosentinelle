// Native marker genes: what the strain's own URA3, HIS3, LEU2, TRP1, LYS2,
// ADE2, MET15, CAN1, ARG4 and LEU5 look like (the list is marker_loci.tsv of the
// database). Reads are aligned to the whole reference genome, not to the loci
// alone, so a read from a paralog or a repeat gets a low mapping quality instead
// of landing in a marker. Only the reads over the loci are kept afterwards,
// plus the depth over the mitochondrial genome and the read count of every
// chromosome, for the mtDNA (rho, growth on glycerol) status, the depth over
// the 2-micron plasmid, whose sequence is added to the reference genome, for
// its status (cir+, cir0) and copy number, and the reads that place the tool
// elements in the genome (ToolPosition).
// The reads are aligned in slices, by MarkerAlignSlice tasks that run in
// parallel, and MarkerAlign merges them. A read is aligned the same way in a
// slice as among all the reads: bowtie2 seeds its pseudo-random choices for
// each read with the read itself.

// bowtie2 indexes of the reference genome, with the 2-micron plasmid, and of
// the tool elements.
process MarkerIndex {
    label 'bowtie2' // see the withLabel: bowtie2 in the nextflow config file
    cache 'true'

    input:
    path reference
    path two_micron_fasta // two_micron.fasta of the database, or NO_2MICRON
    path tools_fasta // tools.fasta of the database, for the positions of the tool elements

    output:
    path "genome_ref.*", emit: genome_index_ch
    path "tools_ref.*", emit: tools_index_ch
    path "*_index_report.txt", emit: marker_index_log_ch

    script:
    """
    touch genome_index_report.txt
    # The 2-micron plasmid is not part of the reference genome: its sequence is
    # added to the genome the reads are aligned to, so that its reads land on it.
    if [[ "${two_micron_fasta}" != "NO_2MICRON" ]]; then
        if [[ "${reference}" == *.gz ]]; then zcat ${reference}; else cat ${reference}; fi | awk 1 > genome.fasta
        awk 1 ${two_micron_fasta} >> genome.fasta
        bowtie2-build --threads ${task.cpus} genome.fasta genome_ref > /dev/null 2>> genome_index_report.txt
        rm genome.fasta
    else
        bowtie2-build --threads ${task.cpus} ${reference} genome_ref > /dev/null 2>> genome_index_report.txt
    fi
    bowtie2-build --threads ${task.cpus} ${tools_fasta} tools_ref > /dev/null 2> tools_index_report.txt
    """
}

// The slices of the reads for MarkerAlignSlice: one per slice_bytes bytes of
// FASTQ (at least one), each with the same number of reads. bowtie2 reads the
// files one after the other, as a single run of reads, and a slice is given by
// the number of reads skipped and the number aligned.
process ReadSlices {
    label 'bash' // see the withLabel: bash in the nextflow config file
    cache 'true'

    input:
    path fastq
    val slice_bytes

    output:
    path "read_slices.txt", emit: read_slices_ch

    script:
    def files = (fastq instanceof List) ? fastq : [fastq]
    """
    set -o pipefail
    for f in ${files.join(' ')} ; do
        if [[ "\$f" == *.gz ]] ; then gzip -dc "\$f" ; else cat "\$f" ; fi | wc -lc
    done | awk -v slice_bytes=${slice_bytes} '
        { reads += int((\$1 + 3) / 4) ; bytes += \$2 }
        END {
            slices = int(bytes / slice_bytes + 0.5)
            if (slices < 1) slices = 1
            size = int((reads + slices - 1) / slices)
            # slice number, reads skipped, reads aligned (0: all the rest, for the last slice)
            for (i = 0 ; i < slices ; i++) print i, i * size, (i < slices - 1 ? size : 0)
        }' > read_slices.txt
    """
}

// The genome alignment of a slice of the reads, and the reads of the slice
// that align on the tool elements.
process MarkerAlignSlice {
    label 'bowtie2' // see the withLabel: bowtie2 in the nextflow config file
    cache 'true'

    input:
    tuple val(slice), val(skip), val(count) // a line of ReadSlices
    path fastq
    path genome_index
    path tools_index

    output:
    path "genome_*.bam", emit: slice_bam_ch
    path "tool_reads_*.tsv", emit: slice_tool_reads_ch
    path "*_align_*.txt", emit: slice_log_ch

    script:
    // Every fastq file, joined for -U: pairing is not needed here, and the
    // second mate is free depth. Mates keep the same read name, which is how
    // ToolPosition finds the mate of a read of a tool element.
    def files = (fastq instanceof List) ? fastq : [fastq]
    def reads = files.join(',')
    def range = "--skip ${skip}" + (count.toLong() > 0 ? " --upto ${count}" : "") // --upto 0 would align no read
    def name = String.format('%05d', slice.toInteger()) // the order of the slices for MarkerAlign
    """
    set -o pipefail # a bowtie2 failure in the pipes below must stop the process
    # Local mode on purpose: reads crossing a junction with foreign sequence (a
    # Ty1, a cassette) are soft-clipped at the junction instead of being dropped,
    # and those clips are the evidence for an insertion.
    bowtie2 --sensitive-local -p ${task.cpus} -x genome_ref -U ${reads} ${range} 2> genome_align_${name}.txt \\
        | samtools sort -@ ${task.cpus} -m 128M -o genome_${name}.bam - # 128 MB at most per thread
    # The reads of the tool elements, aligned end to end on tools.fasta (see MarkerAlign).
    bowtie2 --sensitive -p ${task.cpus} --no-unal --no-hd -x tools_ref -U ${reads} ${range} 2> tools_align_${name}.txt \\
        | awk -F '\\t' -v OFS='\\t' '{print \$1, \$2, \$3, \$4, \$6, \$10}' > tool_reads_${name}.tsv
    """
}

// The slices merged, then what the markers, the mtDNA, the 2-micron plasmid and
// the positions of the tool elements need from the genome alignment.
process MarkerAlign {
    label 'bowtie2' // see the withLabel: bowtie2 in the nextflow config file
    cache 'true'

    input:
    val sample_name
    path fastq
    path slice_bams
    path slice_tool_reads
    path slice_logs
    path index_logs
    path genome_index
    path reference
    path loci_bed
    path mito_genes // mito_genes.tsv of the database, or NO_MITO
    path two_micron_fasta // two_micron.fasta of the database, or NO_2MICRON
    path tools_fasta // tools.fasta of the database, for the positions of the tool elements

    output:
    tuple path("markers.pileup"), path("markers.depth"), path("markers.sam"), emit: marker_evidence_ch
    tuple path("mito.depth"), path("genome.idxstats"), path("read_length.txt"), emit: mito_evidence_ch
    path "two_micron.depth", emit: two_micron_evidence_ch
    tuple path("tool_clips.tsv.gz"), path("tool_reads.tsv"), path("tool_mates.tsv"), path("tool_chunks.tsv"), path("read_length.txt"), path("read_layout.txt"), emit: tool_evidence_ch
    path "marker_align_report.txt", emit: marker_align_log_ch

    script:
    def files = (fastq instanceof List) ? fastq : [fastq]
    def layout = (files.any{ it.name.startsWith('R1__') } ? 'paired' : 'single') +
        (files.any{ it.name.startsWith('S__') } ? ' pieces' : '') // S__: pieces of long sequences
    def bams = ((slice_bams instanceof List) ? slice_bams : [slice_bams]).sort{ it.name }
    def tool_reads = ((slice_tool_reads instanceof List) ? slice_tool_reads : [slice_tool_reads]).sort{ it.name }
    """
    set -o pipefail # a failure in the pipes below must stop the process
    # The alignment summaries of the slices added up, as bowtie2 writes it for
    # unpaired reads, after the warnings.
    summary() {
        awk '
            function percent(part, all) { return all ? 100 * part / all : 0 }
            / reads; of these:\$/ { total += \$1 ; next }
            / were unpaired; of these:\$/ { unpaired += \$1 ; next }
            / aligned 0 times\$/ { none += \$1 ; next }
            / aligned exactly 1 time\$/ { once += \$1 ; next }
            / aligned >1 times\$/ { more += \$1 ; next }
            /^0 reads\$/ || / overall alignment rate\$/ { next }
            { print }
            END {
                printf "%.0f reads; of these:\\n", total
                printf "  %.0f (%.2f%%) were unpaired; of these:\\n", unpaired, percent(unpaired, total)
                printf "    %.0f (%.2f%%) aligned 0 times\\n", none, percent(none, unpaired)
                printf "    %.0f (%.2f%%) aligned exactly 1 time\\n", once, percent(once, unpaired)
                printf "    %.0f (%.2f%%) aligned >1 times\\n", more, percent(more, unpaired)
                printf "%.2f%% overall alignment rate\\n", percent(once + more, total)
            }' "\$@"
    }
    echo -e "\\n\\n<br /><br />\\n\\n### NATIVE MARKER ALIGNMENT for ${sample_name}\\n\\n" > marker_align_report.txt
    cat genome_index_report.txt >> marker_align_report.txt
    summary genome_align_*.txt >> marker_align_report.txt
    two_micron_chrom=""
    if [[ "${two_micron_fasta}" != "NO_2MICRON" ]]; then
        two_micron_chrom=\$(awk 'NR == 1 {print substr(\$1, 2)}' ${two_micron_fasta})
    fi
    if [[ ${bams.size()} -gt 1 ]]; then
        samtools merge -@ ${task.cpus} genome.bam ${bams.join(' ')}
    else
        ln -s ${bams[0]} genome.bam
    fi
    samtools index -@ ${task.cpus} genome.bam
    samtools view -@ ${task.cpus} -b -L ${loci_bed} genome.bam > markers.bam
    samtools index markers.bam

    # mtDNA: reads per chromosome of the reference genome (nuclear depth), the
    # depth of every base of the mitochondrial chromosome, all mapping
    # qualities, and the mean read length.
    samtools idxstats genome.bam | awk -F '\\t' -v plasmid="\${two_micron_chrom}" '\$1 != plasmid' > genome.idxstats
    samtools stats -@ ${task.cpus} genome.bam | awk -F '\\t' '\$1 == "SN" && \$2 == "average length:" {print \$3}' > read_length.txt
    if [[ "${mito_genes}" != "NO_MITO" ]]; then
        mito_chrom=\$(awk -F '\\t' 'NR == 2 {print \$4}' ${mito_genes})
        samtools depth -a -Q 0 -r "\${mito_chrom}" genome.bam > mito.depth
    else
        touch mito.depth
    fi
    # 2-micron plasmid: the depth of every base, read or not (-aa), all mapping
    # qualities (the reads of the two inverted repeats have a low one).
    if [[ -n "\${two_micron_chrom}" ]]; then
        samtools depth -aa -Q 0 -r "\${two_micron_chrom}" genome.bam > two_micron.depth
        echo "reads on the 2-micron plasmid (\${two_micron_chrom}): \$(samtools view -c genome.bam "\${two_micron_chrom}")" >> marker_align_report.txt
    else
        touch two_micron.depth
    fi

    # Tool elements, for their positions in the genome (locate_tools.py, in
    # ToolPosition). (1) The reads of the elements, aligned end to end on
    # tools.fasta (MarkerAlignSlice). (2) In the genome alignment, with a
    # mapping quality of 20 or more: the reads soft-clipped by 15 bases or more
    # (chromosome, junction, side, aligned part, clipped bases, name), the clipped
    # bases being looked up in the elements afterwards, and the reads named like
    # a read of an element (its mate). (3) Pieces of 50 bases of the elements,
    # one every 10 bases, aligned to the genome: where the reference itself
    # carries an element.
    cat tools_index_report.txt >> marker_align_report.txt
    summary tools_align_*.txt >> marker_align_report.txt
    cat ${tool_reads.join(' ')} > tool_reads.tsv
    samtools view -@ ${task.cpus} -F 0x904 -q 20 genome.bam | awk -F '\\t' -v OFS='\\t' -v names=tool_reads.tsv -v min_clip=15 '
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
    path two_micron_depth
    path two_micron_features // two_micron_features.tsv of the database, or NO_2MICRON
    val seq_kind
    val seq_molecule

    output:
    path "marker_genotype.json", emit: marker_json_ch
    path "marker_genotype.tsv", emit: marker_tsv_ch
    path "marker_genotype_report.txt", emit: marker_genotype_log_ch

    script:
    """
    set -o pipefail # a script failure in the pipes to tee below must stop the process
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
    # 2-micron plasmid status (cir+, cir0) and copy number. The nuclear depth
    # leaves out the mitochondrial chromosome, as in mito_status.py.
    two_micron_arg=""
    if [[ "${two_micron_features}" != "NO_2MICRON" ]]; then
        mito_chrom_arg=""
        if [[ "${mito_genes}" != "NO_MITO" ]]; then
            mito_chrom_arg="--mito-chrom \$(awk -F '\\t' 'NR == 2 {print \$4}' ${mito_genes})"
        fi
        two_micron_status.py \\
            --features ${two_micron_features} \\
            --depth ${two_micron_depth} \\
            --idxstats ${idxstats} \\
            --read-length \$(cat ${read_length}) \\
            --seq-kind ${seq_kind} \\
            --seq-molecule ${seq_molecule} \\
            --min-depth ${min_depth} \${mito_chrom_arg} \\
            --out two_micron_status.json |& tee -a marker_genotype_report.txt
        two_micron_arg="--two-micron two_micron_status.json"
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
        --tsv marker_genotype.tsv \${mito_arg} \${two_micron_arg} |& tee -a marker_genotype_report.txt
    """
}
