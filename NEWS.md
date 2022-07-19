# mHMMbayes 0.2.0
## Speed
This release mostly focusses on increasing the speed of the `mHMM()` algorithm. 
* the forward algorithm used in mHMM() is now implemented in c++ using Rcpp to optimize computational speed


## Plots
* ..

## sim_mHMM()
In the function sim_mHMM() used to simulate data for multiple subjects - for which the observation follow a hidden Markov model (HMM) with an multilevel structure - now allows for the simulation of multivariate data.

## Other minor improvements and bug fixes
* fixed error in tutorial vignette, section 'Graphically displaying outcomes'.
* implemented progress bar for mHMM() to indicate progress of the used algorithm


# mHMMbayes 0.1.1

Patch release to solve noLD issuses (tests without long double on x86_64 Linux system) uncoverd by CRAN Package Check Results.

# mHMMbayes 0.1.0 
First (official) version of the package! 
