
process Coverage { // section 24.5 of the labbook 20200707. Warning: USING 5' AND 3' COORDINATES
    label 'bedtools' // see the withLabel: bash in the nextflow config file 
    // publishDir "${out_path}/reports", mode: 'copy', pattern: "cov_report.txt", overwrite: false // inactivated because no cov_report published in "${out_path}/reports" probably because of the parallelization
    //publishDir "${out_path}/files", mode: 'copy', pattern: "*.cov", overwrite: false
    cache 'true'

    input:
    tuple val(ref_name), path(ref), val(file_name), path(bam)

    output:
    path "*_mini.cov", emit: cov_ch // warning: 3 files
    // file "*.cov" // coverage per base if ever required but long process
    path "cov_report.txt", emit: cov_report_ch

    script:
    """
    # bedtools genomecov -d -ibam \${bam} > \${bam.baseName}.cov |& tee cov_report.txt # coverage per base if ever required but long process
    # to add the chr names | awk '{h[\$NF]++}; END { for(k in h) print k, h[k] }' | sort -V > \${bam.baseName}.cov
    bedtools genomecov -bga -ibam ${bam}  > ${ref_name}_mini.cov |& tee cov_report.txt
    # -g \${ref} not required when inputs are bam files
    """
}
