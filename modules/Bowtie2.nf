// Save the config file and the log file for a specific run
process Bowtie2 { // section 24.1 of the labbook 20200707
    label 'bowtie2' // see the withLabel: bash in the nextflow config file 
    cache 'true'

    input:
    val fastq_name
    path fastq
    tuple path(ref), val(ref_name) 

    output:
    tuple val(ref_name), path(ref), val(fastq_name), path("${fastq_name}_vs_${ref_name}.bam"), emit: bowtie2_ch
    path "bowtie2_report.txt", emit: bowtie2_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n###  Bowtie2 indexing of the reference sequence (${ref_name})\\n\\n" > bowtie2_report.txt
    bowtie2-build ${ref} ${ref_name} |& tee -a bowtie2_report.txt
    echo -e "\\n\\n<br /><br />\\n\\n###  Bowtie2 alignment vs ${ref_name}\\n\\n" >> bowtie2_report.txt
    bowtie2 --sensitive -x ${ref_name} -U ${fastq} -t -S ${fastq_name}_vs_${ref_name}.sam |& tee -a tempo.txt
    # --fast: no soft clipping allowed and very sensitive seed alignment
    # -t time displayed
    cat tempo.txt >> bowtie2_report.txt
    sed -i -e ':a;N;\$!ba;s/\\n/\\n<br \\/>/g' tempo.txt
    echo -e "\\n\\n<br /><br />\\n\\n####  samtools conversion\\n\\n" >> bowtie2_report.txt
    samtools view -bh -o tempo.bam ${fastq_name}_vs_${ref_name}.sam |& tee -a bowtie2_report.txt
    samtools sort -o ${fastq_name}_vs_${ref_name}.bam tempo.bam |& tee -a bowtie2_report.txt
    samtools index ${fastq_name}_vs_${ref_name}.bam |& tee -a bowtie2_report.txt
    """
}