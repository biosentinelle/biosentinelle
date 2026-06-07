process Detect_insertions {
    label 'python'
    publishDir path: "${params.out_path}/tsv/insertions", mode: 'copy', pattern: "insertions.json", overwrite: true
    publishDir path: "${params.out_path}/reports", mode: 'copy', pattern: "detect_insertions_report.txt", overwrite: true
    cache 'true'

    input:
    val ref_names
    path ref_bams
    val feat_names
    path feat_bams
    path detect_script

    output:
    path "insertions.json", emit: insertions_json_ch
    path "detect_insertions_report.txt", emit: detect_insertions_report_ch

    script:
    def ref_manifest = ref_names.indices.collect { idx -> "${ref_names[idx]}\t${ref_bams[idx].name}" }.join('\n')
    def feat_manifest = feat_names.indices.collect { idx -> "${feat_names[idx]}\t${feat_bams[idx].name}" }.join('\n')
    """
    #!/bin/bash -ue
    set -o pipefail

    cat > ref_alignments.tsv <<'EOF_REF'
${ref_manifest}
EOF_REF

    cat > feature_alignments.tsv <<'EOF_FEAT'
${feat_manifest}
EOF_FEAT

    python3 ${detect_script} \\
        --refs ref_alignments.tsv \\
        --features feature_alignments.tsv \\
        --output insertions.json \\
        |& tee detect_insertions_report.txt
    """
}
