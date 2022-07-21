# mHMMbayes 0.2.0

## Speed
This release mostly focusses on increasing the speed of the `mHMM()` algorithm. 
* the forward algorithm used in mHMM() is now implemented in c++ using Rcpp to optimize computational speed
* the call to optim() in mHMM() used to create correct scalers for the Metropolis Hasting was computationally very intensive, especially for long sequences of data. In this new version, the log likelihood function of the multinomial distribution is programmed in a more efficient manner. 

## sim_mHMM()
In the function sim_mHMM() used to simulate data for multiple subject - for which the observation follow a hidden Markov model (HMM) with an multilevel structure - now allows for the simulation of multivariate data. The distributions of multiple dependent variables for multivariate data are assumed to be independent given the current hidden state.

## Other minor improvements and bug fixes
* probabilities returned by summary(), obtain_emiss() and obtain_gamma() are now based on the MAP estimates of the intercepts of the multinomial distribution instead of the MAP estimates of the probabilities. This ensures that the returned probabilities sum to 1. 
* implemented progress bar for mHMM() to indicate progress of the used algorithm
* bug fix for plotting subject posterior densities in plot.mHMM(): iterations were set to 10 instead of n_subj. 
* fixed error in tutorial vignette, section 'Graphically displaying outcomes'.


# mHMMbayes 0.1.1
Patch release to solve noLD issuses (tests without long double on x86_64 Linux system) uncoverd by CRAN Package Check Results.

# mHMMbayes 0.1.0 
First (official) version of the package! 
