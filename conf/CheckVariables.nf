workflow CheckVariables {

//////// Checks
//// check of the bin folder
tested_files_bin = [
    "biosentinelle_template.rmd",
    "detect_insertions.py",
    ]
tested_files_bin.each { i1 ->
    if( ! (file("${projectDir}/bin/${i1}").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE ${i1} FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
}
//// end check of the bin folder

//// check of the modules folder
tested_files_modules = [
    "Backup.nf",
    "Bowtie2.nf",
    "CopyLogFile.nf",
    "CountSnps.nf",
    "Coverage.nf",
    "Detect_insertions.nf",
    "Fastq_dump.nf",
    "FindBestRef.nf",
    "Functions.nf",
    "MultiQC.nf",
    "Plot_coverage.nf",
    "Print_report.nf",
    "Print_snp_count.nf",
    "Print_warnings.nf",
    "Prepare_dashboard.nf",
    "Q20.nf",
    "Split_fasta.nf",
    "Unzip.nf",
    "VariantCall.nf",
    "WorkflowParam.nf"
    ]
tested_files_modules.each { i1 ->
    if( ! (file("${projectDir}/modules/${i1}").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE ${i1} FILE MUST BE PRESENT IN THE ./modules FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
}
//// end check of the modules folder

//// check of config file parameters
// Data
if( ! (params.fastq_path in String || params.fastq_path in GString) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID fastq_path PARAMETER IN nextflow.config FILE:\n${params.fastq_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! params.fastq_path =~ /^http.*/){
    // Skip file existence check for remote URLs
    // SRA accessions will be downloaded by fastq-dump, zip files will be processed by Unzip
    if( ! (file(params.fastq_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID fastq_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${params.fastq_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
}
if( ! (params.ref_path in String || params.ref_path in GString) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ref_path PARAMETER IN nextflow.config FILE:\n${params.ref_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else {
    // ref_path can contain multiple space-separated paths
    if( ! (file(params.ref_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ref_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${params.ref_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
}
if( ! (params.features_path in String || params.features_path in GString) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID features_path PARAMETER IN nextflow.config FILE:\n${params.features_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else {
    // features_path can contain multiple space-separated paths
    if( ! (file(params.features_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID features_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${params.features_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
}
//// end check of config file parameters
}
