// Count SNPs from VCF file
process CountSnps {
    label 'bcftools'
    cache 'true'

    input:
    tuple val(ref_name), path(vcf_file)

    output:
    tuple val(ref_name), env(SNP_COUNT), emit: snp_count_ch
    path "snp_count_report.txt", emit: snp_count_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### SNP Counting for ${ref_name}\\n\\n" > snp_count_report.txt
    
    # Count SNPs from VCF
    SNP_COUNT=\$(bcftools view -H -v snps ${vcf_file} 2>/dev/null | wc -l)
    echo "SNPs vs ${ref_name}: \$SNP_COUNT" | tee -a snp_count_report.txt
    """
}
