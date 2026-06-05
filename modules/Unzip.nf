

process Unzip {
    label 'unzip'
    cache 'true'

    input:
    path zip
    val sample_path

    output:
    path "*", emit: unzip_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail # return the exit code of the first nonzero command in the pipeline, not just the last one (important when using tee)
    FILENAME_INI=\$(basename -- "${sample_path}")
    FILE_EXTENSION="\${FILENAME_INI##*.}"
    FILENAME="\${FILENAME_INI%.*}"
    if [[ ! "\${FILE_EXTENSION}" =~ zip ]] ; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nTHE FILE EXTENSION MUST BE \\".zip\\" AND NOT \${FILENAME_INI}\\n\\n========\\n\\n"
        exit 1
    else
        SINGLE_FILE_NAME=\$(unzip -Z1 "${zip}")
        unzip ${zip}
    fi
    rm \$FILENAME_INI # remove the zipped file 
    if [ -d "\$FILENAME" ]; then # if the unzipped file is a directory
        for file in "\$FILENAME"/*.* ; do
            # Check if the file is a regular file and not a directory
            if [ -f "\$file" ] ; then
                # Check if the file does not match the extensions .fasta or .fa
                if [[ ! "\$file" =~ \\.(fasta|fa|fas|fna|txt|seq|faa)\$ ]] ; then
                    echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nALL THE UNZIPPED FILE EXTENSIONS MUST BE fasta, fa, fas, fna, txt, seq OR faa.\\nHERE A FILE IS:\\n\$file.\\n\\n========\\n\\n"
                    exit 1
                fi
            else
                echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nALL THE UNZIPPED FILE EXTENSIONS MUST BE fasta, fa, fas, fna, txt, seq OR faa.\\nHERE A FILE IS:\\n\$file.\\n\\n========\\n\\n"
                exit 1
            fi
        done
        cp "\$FILENAME"/*.* .
        rm -r \$FILENAME
    else
        if [[ ! "\${SINGLE_FILE_NAME}" =~ \\.(fasta|fa|fas|fna|txt|seq|faa)\$ ]] ; then
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nALL THE UNZIPPED FILE EXTENSIONS MUST BE fasta, fa, fas, fna, txt, seq OR faa.\\nHERE A FILE IS:\\n\${SINGLE_FILE_NAME}.\\n\\n========\\n\\n"
            exit 1
        fi
    fi
    """
}