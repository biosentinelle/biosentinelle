// Count the strain panel k-mers in the raw reads and call one allele per site.
// No alignment and no variant calling on the sample: the panel carries the
// k-mer for both alleles, so the reads only have to be streamed once.
process KmerScan {
    label 'python' // see the withLabel: python in the nextflow config file
    publishDir "${out_path}/tsv", mode: 'copy', pattern: "site_calls.tsv", overwrite: false
    cache 'true'

    input:
    val sample_name
    path fastq
    path panel_kmers
    path mating_kmers
    val max_bases
    val min_depth

    output:
    path "site_calls.tsv", emit: site_calls_ch
    path "mating_counts.tsv", emit: mating_counts_ch
    path "kmer_scan_report.txt", emit: kmer_scan_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### PANEL K-MER SCAN for ${sample_name}\\n\\n" > kmer_scan_report.txt

    scan_kmers.py \\
        --kmers ${panel_kmers} \\
        --reads ${fastq} \\
        --out site_calls.tsv \\
        --extra-kmers ${mating_kmers} \\
        --extra-out mating_counts.tsv \\
        --max-bases ${max_bases} \\
        --min-depth ${min_depth} |& tee -a kmer_scan_report.txt
    """
}
