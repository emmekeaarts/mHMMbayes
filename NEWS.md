# mHMM bayes (vary emiss version)
This release focusses on the possibility to model data with different distributions (categorical and continuous) using the multilevel HMM. 
* The function sim_mHMM to simulate multivel HMM data is extended such that data with varying types of emission distributions can be simulated
* new function: mHMM_vary 
     * the structure of the returned object PD_subj is changed. PD_subj contains the subject level paremter estimates and log likelihood over the iterations of the MCMC sampler. Instead of a list containing one matrix per subject containing all parameters, PD_subj is now a list of lists, one list per subject containing the elements `trans_prob`, `cat_emiss`, `cont_emiss`, and `log_likl` which contain the parameter esimates over the iterations of the transition probabilities gamma, the categorical emission distribution, and the continious emission distribuituion, respectiviely, and the log likelihood over the iterations. 


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
