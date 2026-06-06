
process Print_snp_count {
    label 'bash'
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    cache 'true'

    input:
    path tsv

    output:
    path "snp_counts.tsv", emit: final_tsv_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    cp snp_counts.tsv tempo
    """
}
