process MultiQC{
    label "multiqc"
    publishDir "${params.out_path}/reports", mode: 'copy', pattern: "multiqc_report.html", overwrite: false

    input:
    file "*" from fastqc_log_ch.mix(bowtie2_log_ch).mix(krakenreports).collect()

    output:
    file "multiqc_report.html" into multiqc_ch
    file "report.rmd" into log_ch10

    script:
    """
    multiqc . -n multiqc_report.html
    echo -e "\\n\\n<br /><br />\\n\\n###  MultiQC\\n\\n" > report.rmd
    if [[ ${params.system_exec} == "local" ]] ; then
        echo -e "\\n\\nWarning: no Kraken performed when using local run\\n" >> report.rmd
    fi
    echo -e "\\n\\nResults are published in the [Report](./reports/multiqc_report.html) folder\\n" >> report.rmd
    """
}
