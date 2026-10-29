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
    // Reads come from PrepareReads: R1__ and R2__ files are mates, U__ and S__
    // files are unpaired (S__ are pieces of long sequences).
    def files = (fastq instanceof List) ? fastq : [fastq]
    def r1 = files.findAll{ it.name.startsWith('R1__') }.sort{ it.name }
    def r2 = files.findAll{ it.name.startsWith('R2__') }.sort{ it.name }
    def unpaired = files.findAll{ it.name.startsWith('U__') || it.name.startsWith('S__') }
    def reads = (r1 ? "-1 ${r1.join(',')} -2 ${r2.join(',')} " : "") + (unpaired ? "-U ${unpaired.join(',')}" : "")
    """
    echo -e "\\n\\n<br /><br />\\n\\n###  Bowtie2 indexing of the reference sequence (${ref_name})\\n\\n" > bowtie2_report.txt
    bowtie2-build ${ref} ${ref_name} |& tee -a bowtie2_report.txt
    echo -e "\\n\\n<br /><br />\\n\\n###  Bowtie2 alignment vs ${ref_name}\\n\\n" >> bowtie2_report.txt
    set -o pipefail # a bowtie2 failure in the pipe below must stop the process
    bowtie2 --sensitive -x ${ref_name} ${reads} -t 2> tempo.txt | samtools sort -o ${fastq_name}_vs_${ref_name}.bam - 2> samtools_sort.txt
    # --fast: no soft clipping allowed and very sensitive seed alignment
    # -t time displayed
    # sam piped directly into samtools sort, so that no large intermediate .sam or unsorted .bam file is written in the work folder
    cat tempo.txt
    cat tempo.txt >> bowtie2_report.txt
    sed -i -e ':a;N;\$!ba;s/\\n/\\n<br \\/>/g' tempo.txt
    echo -e "\\n\\n<br /><br />\\n\\n####  samtools sort\\n\\n" >> bowtie2_report.txt
    cat samtools_sort.txt >> bowtie2_report.txt
    echo -e "\\n\\n<br /><br />\\n\\n####  samtools index\\n\\n" >> bowtie2_report.txt
    samtools index ${fastq_name}_vs_${ref_name}.bam |& tee -a bowtie2_report.txt
    """
}