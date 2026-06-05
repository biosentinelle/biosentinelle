
process Print_warnings {
    label 'bash'

    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "warnings.txt", overwrite: false
    cache 'true'

    input:
    path warning_ch

    output:
    path "warnings.txt", emit: final_warning_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    cp warnings_collect.txt warnings.txt
    """
}
