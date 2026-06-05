
process CopyLogFile{
    label 'bash'
    cache 'true'

    input:
    path log_file_name
    val out_path 

    script:
    """
    cp "${log_file_name}" "${out_path}/reports/${log_file_name.getName()}"
    """
}


