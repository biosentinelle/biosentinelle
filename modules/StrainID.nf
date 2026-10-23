// Combine the panel allele calls and the tool element scores into the strain
// identity result, then render the interactive HTML report.
process StrainScore {
    label 'python' // see the withLabel: python in the nextflow config file
    publishDir "${out_path}/tsv", mode: 'copy', pattern: "strain_id.json", overwrite: false
    cache 'true'

    input:
    val sample_name
    path site_calls
    path panel_sites
    path tools_manifest
    path tool_coverage
    path marker_json
    val organism
    val reference_name
    val error_rate
    val min_concordance

    output:
    path "strain_id.json", emit: strain_json_ch
    path "strain_score_report.txt", emit: strain_score_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### STRAIN SCORING for ${sample_name}\\n\\n" > strain_score_report.txt

    score_strains.py \\
        --sites ${panel_sites} \\
        --calls ${site_calls} \\
        --tools ${tools_manifest} \\
        --tool-coverage ${tool_coverage} \\
        --markers ${marker_json} \\
        --sample "${sample_name}" \\
        --organism "${organism}" \\
        --reference "${reference_name}" \\
        --error-rate ${error_rate} \\
        --min-concordance ${min_concordance} \\
        --out strain_id.json |& tee -a strain_score_report.txt
    """
}

process StrainReport {
    label 'python' // see the withLabel: python in the nextflow config file
    publishDir "${out_path}", mode: 'copy', pattern: "strain_identity_report.html", overwrite: false
    cache 'true'

    input:
    val sample_name
    path strain_json

    output:
    path "strain_identity_report.html", emit: strain_report_ch

    script:
    """
    render_report.py \\
        --json ${strain_json} \\
        --out strain_identity_report.html \\
        --title "BioSentinelle strain identity - ${sample_name}"
    """
}
