# Configuration Files

These files contain variables for repository path, output data path, data paths, etc. 
The variables must be changed to be computer/project specific. 

Specifically, the following changes are required within the [config.R](config.R) file:

1. `PROJECT_PATH` (location of the `microhomology` repository) and `OUTPUT_PATH` (location where output files should be stored) variables
2. Required data file path (`TCR_REPERTOIRE_DATA_igor` variable for the IGoR processed training data set)
3. Path to the full germline gene names and sequences from IMGT (`WHOLE_NUCSEQS_igor` variable)
