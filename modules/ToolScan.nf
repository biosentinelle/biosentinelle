// Score each element of the tool database by breadth of coverage and identity.
// Alignment is done by the existing Bowtie2 process, so these two steps start
// from its bam and stay inside the containers that already exist.

// Extract per-base depth and per-read mismatch counts from the alignment.
process ToolDepth {
    label 'samtools' // see the withLabel: samtools in the nextflow config file
    cache 'true'

    input:
    tuple val(ref_name), path(ref), val(file_name), path(bam)

    output:
    tuple val(ref_name), path("depth.txt"), path("alignments.txt"), emit: tool_depth_ch
    path "tool_depth_report.txt", emit: tool_depth_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### TOOL DEPTH for ${ref_name}\\n\\n" > tool_depth_report.txt

    samtools sort -o sorted.bam ${bam} |& tee -a tool_depth_report.txt
    samtools index sorted.bam

    # -a keeps the zero-depth positions. Without it, breadth would be computed
    # over the covered part only, and every element would look fully covered.
    samtools depth -a sorted.bam > depth.txt

    # Reference, CIGAR and NM per aligned read. NM carries the edit distance,
    # which is what turns breadth into an identity measurement.
    samtools view -F 4 sorted.bam | awk -F'\\t' '{
        nm = "NA"
        for (i = 12; i <= NF; i++) if (\$i ~ /^NM:i:/) { split(\$i, a, ":"); nm = a[3] }
        print \$3 "\\t" \$6 "\\t" nm
    }' > alignments.txt

    echo "aligned records: \$(wc -l < alignments.txt)" >> tool_depth_report.txt
    """
}

// Combine depth and identity into one score per element.
process ToolScore {
    label 'python' // see the withLabel: python in the nextflow config file
    publishDir "${out_path}/tsv", mode: 'copy', pattern: "tool_coverage.tsv", overwrite: false
    cache 'true'

    input:
    tuple val(ref_name), path(depth), path(alignments)
    path tools_manifest

    output:
    path "tool_coverage.tsv", emit: tool_coverage_ch
    path "tool_score_report.txt", emit: tool_score_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### TOOL ELEMENT SCORING\\n\\n" > tool_score_report.txt

    tool_coverage.py \\
        --depth ${depth} \\
        --alignments ${alignments} \\
        --manifest ${tools_manifest} \\
        --out tool_coverage.tsv |& tee -a tool_score_report.txt
    """
}
