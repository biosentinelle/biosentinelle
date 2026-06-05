workflow CheckVariables {

//////// Checks
//// check of the bin folder
tested_files_bin = [
    "biosentinelle_template.rmd", 
    ]
for(i1 in tested_files_bin){
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
    "Fastq_dump.nf", 
    "Functions.nf",  
    "Print_report.nf", 
    "Print_warnings.nf", 
    "Unzip.nf", 
    "WorkflowParam.nf"
    ]
for(i1 in tested_files_modules){
    if( ! (file("${projectDir}/modules/${i1}").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE ${i1} FILE MUST BE PRESENT IN THE ./modules FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
}
//// end check of the modules folder

//// check of config file parameters
// Data
if( ! (fastq_path in String || fastq_path in GString) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID fastq_path PARAMETER IN nextflow.config FILE:\n${fastq_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (file(fastq_path).exists()) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID fastq_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${fastq_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
}
if( ! (ref_path in String || ref_path in GString) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ref_path PARAMETER IN nextflow.config FILE:\n${ref_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (file(ref_path).exists()) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ref_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${ref_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
}
//// end check of config file parameters
}
