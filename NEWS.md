# mHMMbayes 0.2.0

# mHMMbayes (development version)
## Speed
This release mostly focusses on increasing the speed of the `mHMM()` algorithm. 
* the forward algorithm used in mHMM() is now implemented in c++ using Rcpp to optimize computational speed
* the call to optim() in mHMM() used to create correct scalers for the Metropolis Hasting was computationally very intensive, especially for long sequences of data. In this new version, the output returned by optim() is saved, and optim() is only called when the required output was not generated yet. 

## Plots
* ..

## sim_mHMM()
In the function sim_mHMM() used to simulate data for multiple subject - for which the observation follow a hidden Markov model (HMM) with an multilevel structure - the following updates have been made: 

* sim_mHMM() now allows for the simulation of continuous data. That is, the emission distribution follows a normal (i.e., Gaussian) distribution.
* sim_mHMM() now allows for the simulation of multivariate data.
* The distributions of multiple dependent variables for multivariate data are assumed to be independent, and all distributions for one dataset have to be of the same type (either categorical or continuous).  

## Other minor improvements and bug fixes
* fixed error in tutorial vignette, section 'Graphically displaying outcomes'.
* implemented progress bar for mHMM() to indicate progress of the used algorithm


# mHMMbayes 0.1.1
Patch release to solve noLD issuses (tests without long double on x86_64 Linux system) uncoverd by CRAN Package Check Results.

# mHMMbayes 0.1.0 
First (official) version of the package! 
