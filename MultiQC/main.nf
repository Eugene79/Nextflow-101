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

// Create FASTP process
process FASTP {
    container 'biocontainers/fastp:v0.20.1_cv1'
    
    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path("*_trimmed.fastq.gz"), emit: trimmed_reads
    path "*.json", emit: fastp_json
    path "*.html", emit: fastp_html

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
        -o ${name}_R1_trimmed.fastq.gz -O ${name}_R2_trimmed.fastq.gz \
        -j ${name}_fastp.json -h ${name}_fastp.html
    """
}

// Define the main workflow
workflow {
    // Create a channel for input reads
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    read_pairs_ch.view()

    
    // Run FastQC on raw reads
    fastqc_results = FASTQC(read_pairs_ch)

    // Run FASTP on raw reads for filtering and trimming
    fastp_results = FASTP(read_pairs_ch)
}