process Prepare_dashboard {
    label 'bash'
    publishDir path: "${out_path}", mode: 'copy', pattern: "dashboard", overwrite: true
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "dashboard_report.txt", overwrite: true
    cache 'true'

    input:
    path insertions_json
    path dashboard_assets

    output:
    path "dashboard", emit: dashboard_dir_ch
    path "dashboard_report.txt", emit: dashboard_report_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail

    mkdir -p dashboard
    cp ${dashboard_assets} dashboard/
    cp ${insertions_json} dashboard/insertions.json
    printf 'window.BIOSENTINEL_INSERTIONS = ' > dashboard/insertions-data.js
    cat ${insertions_json} >> dashboard/insertions-data.js
    printf ';\\n' >> dashboard/insertions-data.js

    echo "Dashboard prepared in ${out_path}/dashboard" > dashboard_report.txt
    """
}
