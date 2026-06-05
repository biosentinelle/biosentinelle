// Save the config file and the log file for a specific run
process Bowtie2 { // section 24.1 of the labbook 20200707
    label 'bowtie2' // see the withLabel: bash in the nextflow config file 
    cache 'true'

    input:
    val file_name
    val ref_name
    path fastq
    path ref

    output:
    path "${file_name}_bowtie2.bam", emit: bowtie2_ch
    path "bowtie2_report.txt", emit: bowtie2_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n###  Bowtie2 indexing of the reference sequence\\n\\n" >> bowtie2_report.txt
    bowtie2-build ${ref} ${ref_name} |& tee -a bowtie2_report.txt
    echo -e "\\n\\n<br /><br />\\n\\n###  Bowtie2 alignment\\n\\n" > report.rmd
    echo -e "\\n\\n<br /><br />\\n\\n###  Bowtie2 alignment\\n\\n" >> bowtie2_report.txt
    bowtie2 --very-sensitive -x ${ref_name} -U ${fastq} -t -S ${file_name}_bowtie2.sam |& tee -a tempo.txt
    # --very-sensitive: no soft clipping allowed and very sensitive seed alignment
    # -t time displayed
    cat tempo.txt >> bowtie2_report.txt
    sed -i -e ':a;N;\$!ba;s/\\n/\\n<br \\/>/g' tempo.txt
    cat tempo.txt >> report.rmd
    echo -e "\\n\\n<br /><br />\\n\\n####  samtools conversion\\n\\n" >> bowtie2_report.txt
    # samtools faidx ${ref}
    samtools view -bh -o tempo.bam ${file_name}_bowtie2.sam |& tee -a bowtie2_report.txt
    samtools sort -o ${file_name}_bowtie2.bam tempo.bam |& tee -a bowtie2_report.txt
    samtools index ${file_name}_bowtie2.bam |& tee -a bowtie2_report.txt
    """
}