// Count SNPs from VCF file
process SexCount {
    label 'r_ig_clustering'
    cache 'true'

    input:
    tuple val(ref_name), path(vcf_file)

    output:
    tuple val(ref_name), path("${ref_name}.tsv"), emit: sex_count_ch
    path "sex_count_report.txt", emit: sex_count_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### sex geno Counting for ${ref_name}\\n\\n" > sex_count_report.txt
    Rscript -e '
        df <- read.table(${vcf_file}, header = TRUE, stringsAsFactors = FALSE, sep = "\\t") # does not take the header
        df[, 6] <- factor(df[, 6], levels = c("0/0", "0/1", "1/0", "1/1"))
        tab <- table(df[, 6])
        df2 <- data.frame("${ref_name}", tab)
        names(df2) <- c("SEX, "GENO", "FREQ")
        write.table(df2, file = paste0(${ref_name}, ".tsv"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
' | tee -a sex_count_report.txt
    """
}
