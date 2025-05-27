// Default parameter input
params.str = "Hello world!"

// splitString process
process splitString {
    publishDir "results/lower"
    
    input:
    val x
    
    output:
    path 'chunk_*'

    script:
    """
    printf '${x}' | split -b 6 - chunk_
    """
}

// convertToUpper process
process convertToUpper {
    publishDir "results/upper"
    tag "$y"

    input:
    path y

    output:
    path 'upper_*'

    script:
    """
    cat $y | tr '[a-z]' '[A-Z]' > upper_${y}
    """
}

// Workflow block
workflow {
    ch_str = channel.of(params.str)     // Create a channel using parameter input
    ch_chunks = splitString(ch_str)     // Split string into chunks and create a named channel
    convertToUpper(ch_chunks.flatten()) // Convert lowercase letters to uppercase letters
}