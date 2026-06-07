process Sex2tsv {
    label 'bcftools'
    cache 'true'

    input:
    tuple val(ref_name), path(vcf_file)

    output:
    tuple val(ref_name), path("sex_geno.tsv"), emit: var_count_ch
    path "var_count_report.txt", emit: var_count_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### VAR Counting for ${ref_name}\\n\\n" > var_count_report.txt
    bcftools view -H ${vcf_file} 2>/dev/null > sex_geno.tsv | tee -a var_count_report.txt
    """
}