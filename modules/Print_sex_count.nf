
process Print_sex_count {
    label 'bash'
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    cache 'true'

    input:
    path tsv

    output:
    path "sex.tsv", emit: final_sex_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    cp sex.tsv tempo
    """
}
