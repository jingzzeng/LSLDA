This folder contains 10 files:

- lslda.R: The main file. It contains the main function which implements the LSLDA algorithm. The validation functions for tuning parameter selection are also included.
- utility.R: Contains some utility functions.
- simulation1.R-simulation7.R: These files implements LSLDA on simulated data generated from models (M1)-(M7) (Section 5), where p = 3000. 
- README.txt: This file.

Note:
1. The comments in files "simulation1.R" to "simulation7.R" are similar, we refer the users to "simulation1.R" for the helpful comments.
2. The files "simulation1.R" to "simulation7.R" only run 2 replicates to save the execution time. The users can specify the variable "times" in the codes to change the number of replicates.