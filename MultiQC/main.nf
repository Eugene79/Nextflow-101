nextflow.enable.dsl = 2

params.reads = "${params.input_fastqs}/*{_R1,_R2}*.fastq.gz"

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
    container 'quay.io/biocontainers/fastp:v0.20.1_cv1'
    
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

// MultiQC process
process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0'

    publishDir 'out', mode: 'copy'
    
    input:
    path '*'

    output:
    path "multiqc_report.html", emit: report, optional: true
    path "multiqc_data", optional: true

    script:
    """
    multiqc .
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

    // Collect all QC reports
    multiqc_files = fastqc_results.fastqc_results
        .mix(fastp_results.fastp_json)
        .mix(fastp_results.fastp_html)
        .collect()

    multiqc_files.view()

    // Run MultiQC
    MULTIQC(multiqc_files)

}
