// Save the config file and the log file for a specific run
process Fastq_dump { // section 24.1 of the labbook 20200707
    label 'fastq_dump' // see the withLabel: bash in the nextflow config file 
    cache 'true'

    input:
    val fastq_path

    output:
    path "*.fastq", emit: fastq_dump_ch
    path "fastq_dump_report.txt", emit: fastq_dump_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### SRA decompression\\n\\n" >> fastq_dump_report.txt
    fastq-dump --split-3 ${fastq_path} |& tee -a fastq_dump_report.txt
    """
}