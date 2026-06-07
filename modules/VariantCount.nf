// Count variants from VCF file
process VariantCount {
    label 'bcftools'
    cache 'true'

    input:
    tuple val(ref_name), path(vcf_file)

    output:
    tuple val(ref_name), env(VAR_COUNT), emit: var_count_ch
    path "var_count_report.txt", emit: var_count_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### VAR Counting for ${ref_name}\\n\\n" > var_count_report.txt
    
    # Count VARs from VCF
    VAR_COUNT=\$((\$(bcftools view -H ${vcf_file} 2>/dev/null | wc -l) - 1))
    echo "VARs vs ${ref_name}: \$VAR_COUNT" | tee -a var_count_report.txt
    """
}
