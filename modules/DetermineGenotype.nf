// Determine genotype (MATA or MATALPHA) based on SNP counts
process DetermineGenotype {
    label 'bash'
    cache 'true'

    input:
    tuple val(ref_name), val(snp_count)

    output:
    path "genotype_report.txt", emit: genotype_report_ch

    script:
    """
    echo "Genotype: ${ref_name}" > genotype_report.txt
    echo "SNP Count: ${snp_count}" >> genotype_report.txt
    echo "" >> genotype_report.txt
    
    # Map the reference name to genotype
    if [[ "${ref_name}" == "MATA" ]]; then
        echo "DETECTED GENOTYPE: MATA" >> genotype_report.txt
    elif [[ "${ref_name}" == "MATALPHA" ]]; then
        echo "DETECTED GENOTYPE: MATALPHA" >> genotype_report.txt
    else
        echo "DETECTED GENOTYPE: ${ref_name}" >> genotype_report.txt
    fi
    """
}
