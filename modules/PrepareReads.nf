// Input normalization: whatever the sample format, the rest of the pipeline
// receives FASTQ files whose name prefix says how to align them (see
// bin/prepare_reads.py).

// BAM, SAM and CRAM back to FASTQ. The original alignment is not reused: it is
// against a reference of somebody else's choosing, and the pipeline aligns
// against its own. Secondary and supplementary records are dropped (samtools
// fastq default), so each read comes out once.
process AlignmentToFastq {
    label 'samtools' // see the withLabel: samtools in the nextflow config file
    cache 'true'

    input:
    path aln
    path cram_ref // a placeholder file named NO_CRAM_REF when none is given

    output:
    path "*.fastq.gz", emit: aln_fastq_ch
    path "aln_to_fastq_report.txt", emit: aln_to_fastq_log_ch

    script:
    def stem = aln.name.replaceAll(/\.(bam|sam|cram)$/, '')
    // Without a reference, samtools looks the CRAM reference up by its MD5
    // (REF_PATH, or the EBI server, which needs network access).
    def ref_opt = cram_ref.name == 'NO_CRAM_REF' ? '' : "--reference ${cram_ref}"
    """
    #!/bin/bash -ue
    set -o pipefail
    echo -e "\\n\\n<br /><br />\\n\\n### ${aln.name} -> FASTQ\\n\\n" > aln_to_fastq_report.txt
    # collate puts mates next to each other, which samtools fastq needs to
    # split pairs. Unpaired reads (long reads, single end) go to the _0 file.
    samtools collate ${ref_opt} -O -u ${aln} tmp_collate \\
        | samtools fastq -n -1 ${stem}_1.fastq.gz -2 ${stem}_2.fastq.gz -0 ${stem}_0.fastq.gz -s ${stem}_single.fastq.gz - \\
        |& tee -a aln_to_fastq_report.txt
    # Remove the outputs that stayed empty.
    for f in ${stem}_1.fastq.gz ${stem}_2.fastq.gz ${stem}_0.fastq.gz ${stem}_single.fastq.gz ; do
        if [[ -f \$f && -z \$(gzip -cd \$f | head -c 1 || true) ]] ; then rm \$f ; fi
    done
    if ! ls *.fastq.gz > /dev/null 2>&1 ; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nNO READ IN ${aln.name}\\n\\n========\\n\\n"
        exit 1
    fi
    """
}

// Any sequence file to the align/ and scan_only/ read sets.
process PrepareReads {
    label 'python' // see the withLabel: python in the nextflow config file
    publishDir "${out_path}/tsv", mode: 'copy', pattern: "prepare_reads.tsv", overwrite: false
    cache 'true'

    input:
    val sample_name
    path seq_files
    val seq_kind
    val max_length
    val frag_length

    output:
    path "align/*", emit: align_ch
    path "scan_only/*", optional: true, emit: scan_only_ch
    path "prepare_reads.tsv", emit: prepare_tsv_ch
    path "prepare_reads_report.txt", emit: prepare_log_ch

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n### INPUT SEQUENCES of ${sample_name}\\n\\n" > prepare_reads_report.txt
    prepare_reads.py ${seq_files} \\
        --kind ${seq_kind} \\
        --max-length ${max_length} \\
        --frag-length ${frag_length} \\
        --report prepare_reads.tsv |& tee -a prepare_reads_report.txt
    cat prepare_reads.tsv >> prepare_reads_report.txt
    """
}

// Reference, feature and mating type files: GenBank and EMBL to FASTA, FASTA
// passed through. All the files of one folder in one task.
process ToFasta {
    label 'python' // see the withLabel: python in the nextflow config file
    cache 'true'

    input:
    path ref_files

    output:
    path "fasta/*.fasta", emit: fasta_ch

    script:
    """
    prepare_reads.py --to-fasta --out-dir fasta ${ref_files}
    """
}
