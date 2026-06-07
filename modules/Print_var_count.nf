
process Print_var_count {
    label 'bash'
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    cache 'true'

    input:
    path tsv

    output:
    path "variant_counts.tsv", emit: final_tsv_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    cp variant_counts.tsv tempo
    """
}
