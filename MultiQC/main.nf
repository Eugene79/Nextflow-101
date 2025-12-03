// Pipeline parameters (example)
params.input_file = "input.txt"

// Processes (example)
process STEP1 {
  input:
    path data

  output:
    path "output.txt"

  script:
    """
    echo "$data" > output.txt
    """
}

// Workflow (example)
workflow {
  //create an input channel
  ch_input = channel.fromPath(params.input_file)  
  // Connect the input channel to the process
  STEP1( ch_input )
}