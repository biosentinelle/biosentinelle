// Download of a run accession of SRA, ENA or DDBJ
process Fastq_dump { // section 24.1 of the labbook 20200707
    label 'fastq_dump' // see the withLabel: fastq_dump in the nextflow config file
    cache 'true'

    input:
    val fastq_path

    output:
    path "*.fastq", emit: fastq_dump_ch
    path "fastq_dump_report.txt", emit: fastq_dump_log_ch

    script:
    """
    set -o pipefail # a failed download stops the run, instead of passing on truncated files
    echo -e "\\n\\n<br /><br />\\n\\n### SRA decompression\\n\\n" >> fastq_dump_report.txt
    fastq-dump --split-3 ${fastq_path} |& tee -a fastq_dump_report.txt
    """
}

// Gzipped copy of each downloaded file in dataset/<accession>/, where main.nf
// is (folder ignored by git). This folder can then be given as fastq_path, to
// analyze the run again without downloading it. The analysis does not wait for
// this copy: it uses the uncompressed files of Fastq_dump.
process Fastq_save {
    label 'bash' // see the withLabel: bash in the nextflow config file
    publishDir "${projectDir}/dataset/${accession}", mode: 'copy', overwrite: false
    cache 'true'

    input:
    val accession
    path fastq

    output:
    path "${fastq}.gz"

    script:
    """
    gzip -c ${fastq} > ${fastq}.gz
    """
}
