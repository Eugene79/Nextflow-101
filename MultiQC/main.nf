nextflow.enable.dsl = 2

// Define container directives
process FASTQC {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    
    input:
    tuple val(name), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

// Define the main workflow
workflow {
    // Create a channel for input reads
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    read_pairs_ch.view()

    
    // Run FastQC on raw reads
    fastqc_results = FASTQC(read_pairs_ch)

}