
process Q20 { 
    label 'samtools' // see the withLabel: bash in the nextflow config file  
    cache 'true'

    input:
    tuple val(ref_name), path(ref), val(file_name), path(bam)

    output:
    tuple val(ref_name), path(ref), val(file_name), path("${file_name}_q20_dup.bam"), emit: q20_ch
    val "read_nb_before", emit: q20_bef_read_nb_ch
    val "read_nb_after", emit: q20_after_read_nb_ch
    path "q20_report.txt", emit: q20_report_ch

    script:
    """
    samtools view -q 20 -b ${bam} > ${file_name}_q20_dup.bam |& tee q20_report.txt
    samtools index ${file_name}_q20_dup.bam
    echo -e "\\n\\n<br /><br />\\n\\n###  Q20 filtering\\n\\n" > report.rmd
    read_nb_before=\$(samtools view ${bam} | wc -l | cut -f1 -d' ') # -h to add the header
    read_nb_after=\$(samtools view ${file_name}_q20_dup.bam | wc -l | cut -f1 -d' ') # -h to add the header
    """
}