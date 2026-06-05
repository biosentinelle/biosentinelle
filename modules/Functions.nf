// modules/subscribe_helpers.nf
def reportEmptyProcess(process_name, ch) {
    ch.count().subscribe { n ->
        if (n == 0)
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING ${process_name}\n\n========\n\n"
    }
}

def copyLogFile(logfilename, ch, out_path) {
    ch.collectFile(name: "${logfilename}").subscribe { f ->
        f.copyTo("${out_path}/reports") }
}