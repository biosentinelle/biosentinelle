// Variant calling using bcftools mpileup and bcftools call
process VariantCall {
    label 'bcftools'
    cache 'true'

    input:
    tuple val(ref_name), path(ref_file), val(bam_name), path(bam_file)

    output:
    tuple val(ref_name), path("${ref_name}_variants.vcf.gz"), emit: vcf_ch
    path "variant_call_report.txt", emit: variant_call_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### Variant Calling vs ${ref_name}\\n\\n" > variant_call_report.txt
    
    # Index the BAM file if not already indexed
    # samtools index ${bam_file} 2>&1 | tee -a variant_call_report.txt
    
    # Call variants using bcftools
    echo -e "\\n### bcftools mpileup and call\\n" >> variant_call_report.txt
    bcftools mpileup -f ${ref_file} -q 20 -Q 20 -Ou ${bam_file} | \
    bcftools call -mv -Oz -o ${ref_name}_variants.vcf.gz 2>&1 | tee -a variant_call_report.txt
    
    # Index the VCF
    bcftools index ${ref_name}_variants.vcf.gz 2>&1 | tee -a variant_call_report.txt
    
    echo -e "\\nVariant calling complete for ${ref_name}" >> variant_call_report.txt
    """
}
