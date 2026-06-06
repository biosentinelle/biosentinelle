
process Split_fasta {
    label 'bash'
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{split_fasta.log}", overwrite: false

    cache 'true'

    input:
    path fs_ch // no parallelization

    output:
    path "*.fasta", emit: split_fasta_ch
    path "split_fasta.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    # protect linked source
    FILENAME=\$(basename -- ${fs_ch}) # recover a file name without path
    cp -Lr ${fs_ch} "./TEMPO.caca" # to have the hard file, not the symlink, because modifications will be performed inside
    chmod 777 TEMPO.caca
    rm \$FILENAME # remove the initial file to avoid to send it into the channel
    cp -rp TEMPO.caca "\$FILENAME" # -p for preserve permissions
    rm TEMPO.caca
    # end protect linked source

    # variables
    FILE=\${FILENAME%.*} # file name without extension
    MODIF_FILE="\${FILE// /_}" # spaces in the name file replaced by _
    FILE="\${MODIF_FILE}"
    FILE_EXTENSION="\${FILENAME##*.}" #  ## means "delete the longest regex starting at the beginning of the tested string". If nothing, delete nothing. Thus ##*. means delete the longest string finishing by a dot. Use # instead of ## for "delete the shortest regex starting at the beginning of the tested string"
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a split_fasta.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a split_fasta.log
    # end variables

    # checks
    if [[ ! "\${FILE_EXTENSION}" =~ fasta|fa|fas|fna|txt|seq ]] ; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nINVALID FILE EXTENSION IN THE sample_path PARAMETER OF THE nextflow.config FILE:\\n\${FILENAME}\\nMUST BE fasta|fa|fas|fna|txt|seq|faa\\n\\n========\\n\\n"
        exit 1
    fi
    # end checks
    # remove carriage returns
    sed 's/\\r\$//' \${FILENAME} > tempo_file
    # end remove carriage returns
    # remove \\n in the middle of the sequence 
    awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){s=substr(\$0,1,1); rest=substr(\$0,2); gsub(/[^a-zA-Z0-9]/,"_",rest) ; if(NR>1){print "\\n"} ; if (length(rest) > 100){print s substr(rest,1,100)"\\n"}else{print s rest"\\n"}} else {print \$0 ; next}}END{print "\\n"}' tempo_file > \${FILE}.ttt
    # gsub(/[^a-zA-Z0-9]/,"_",rest) replace any weird chars in the first line by a single underscore
    # \${FILE}.ttt is a trick to do not use .fa and modify the initial file due to the link in the nextflow work folder
    # end remove \\n in the middle of the sequence 

    TEMPO=\$(wc -l \${FILE}.ttt | cut -f1 -d' ')
    if read -n 1 char <"\${FILE}.ttt"; [[ \$char != ">" || \$TEMPO != 2 ]]; then
        awk '/^>/ {F=substr(\$1,2)".fasta"} {print >> F}' "\${FILE}.ttt"
    else
        cp "\${FILE}.ttt" "\${FILE}.fasta"
    fi
    """
}