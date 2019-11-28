
# mHMMbayes (development version)
## Speed
This release mostly focusses on increasing the speed of the `mHMM()` algorithm. 
* the forward algorithm used in mHMM() is now implemented in c++ using Rcpp to optimize computational speed
* the call to optim() in mHMM() used to create correct scalers for the Metropolis Hasting was computationally very intensive, especially for long sequences of data. In this new version, the output returned by optim() is saved, and optim() is only called when the required output was not generated yet. 

## Plots
* ..

## Other minor improvements and bug fixes
* fixed error in tutorial vignette, section 'Graphically displaying outcomes'.
* implemented progress bar for mHMM() to indicate progress of the used algorithm


# mHMMbayes 0.1.1
Patch release to solve noLD issuses (tests without long double on x86_64 Linux system) uncoverd by CRAN Package Check Results.

# mHMMbayes 0.1.0 
First (official) version of the package! 
