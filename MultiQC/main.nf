// Pipeline parameters (example)
params.input_file = "input.txt"

// Channels (example)
channel {
  fromFile(params.input_file)
}

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
  // Connect the input channel to the process
  STEP1( ch_input )
}