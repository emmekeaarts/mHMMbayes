# mHMMbayes 0.2.0

## Speed
A major improvement in this release is the increased speed of the `mHMM()` algorithm. 
* the forward algorithm used in mHMM() is now implemented in c++ using Rcpp to optimize computational speed
* the call to optim() in mHMM() used to create correct scalers for the Metropolis Hasting was computationally very intensive, especially for long sequences of data. In this new version, the log likelihood function of the Multinomial distribution is programmed in a more efficient manner, and obtaining the Hessian based on the outcomes of optim() is done more efficiently.  

## Manually specifying hyper-prior distribution paramter values
Two new functions to manually specify hyper-prior distribution parameter values for the multilevel hidden Markov model are introduced:

* prior_emiss_cat(): for manually specifying hyper-prior distribution parameter values for the categorical emission distribution(s), creating an object of class 'mHMM_prior_emiss'. 
* prior_gamma(): for manually specifying hyper-prior distribution parameter values for the transition probability matrix gamma, creating an object of class 'mHMM_prior_gamma'. 

Using manually specified hyper-prior distribution parameter values in the function mHMM() is as of now thus done by inputting an object of the class 'mHMM_prior_emiss' and/or 'mHMM_prior_gamma' for the input parameters emiss_hyp_prior and gamma_hyp_prior, respectively, created by the above functions. Note that manually specifying hyper-prior distribution parameter values is optional, default values are available for all parameters. 

## Transforming a set of probabilities to Multinomial logit regression intercepts and vice versa
Manually specifying hyper-prior distribution parameter values is done on the logit domain. That is, the hyper-priors are on the intercepts (and, if subject level covariates are used, regression coefficients) of the Multinomial logit model used to accommodate the multilevel framework of the data, instead of on the probabilities directly. As logit domain might be more unfamiliar to the user compared to the probability domain, two functions are introduced to aid the user: 

* prob_to_int(): transforms a set of state transition or categorical emission observation probabilities to the corresponding Multinomial logit regression intercepts.
* int_to_prob(): transforms a set of Multinomial logit regression intercepts to the corresponding state transition or categorical emission observation probabilities 

## Manually specifying settings of the proposal distribution of the Random Walk Metropolis sampler 
Two new functions to manually specify settings of the proposal distribution of the Random Walk (RW) Metropolis sampler for the multilevel hidden Markov model are introduced:

* pd_RW_emiss_cat(): for manually specifying setting of the RW proposal distribution for the categorical emission distribution(s), creating an object of class 'mHMM_pdRW_emiss'. 
* pd_RW_gamma(): for manually specifying setting of the RW proposal distribution for the transition probability matrix gamma, creating an object of class 'mHMM_pdRW_gamma'.  

Using manually specified settings of the proposal distribution of the Random Walk (RW) Metropolis sampler in the function mHMM() is as of now thus done by inputting an object of the class 'mHMM_pdRW_emiss' and/or 'mHMM_pdRW_gamma' for the input parameters emiss_sampler and gamma_sampler, respectively, created by the above functions. Note that manually specifying setting of the RW proposal distribution is optional, default values are available for all parameters. 

## sim_mHMM()
In the function sim_mHMM() used to simulate data for multiple subject - for which the observation follow a hidden Markov model (HMM) with an multilevel structure 

* now allows for the simulation of multivariate data. The distributions of multiple dependent variables for multivariate data are assumed to be independent given the current hidden state.
* as such, the input parameter 'gen' is introduced, similar as used in e.g., the function mHMM(). gen contains the elements m, n_dep and q_emiss. gen replaces the input parameters m and q_emiss. In this version, however, using m and q_emiss and not specifying gen or n_dep issues a warning and is thus still allowed. 

## Other minor (quite technical) improvements and bug fixes
* implemented progress bar for mHMM() to indicate progress of the used algorithm
* probabilities returned by summary(), obtain_emiss() and obtain_gamma() are now based on the MAP estimates of the intercepts of the Multinomial distribution instead of the MAP estimates of the probabilities. This ensures that the returned probabilities sum to 1. 
* bug fix for plotting subject posterior densities in plot.mHMM(): iterations were set to 10 instead of n_subj. 
* fixed error in tutorial vignette, section 'Graphically displaying outcomes'.
* for subject specific transition probability matrix gamma and emission distributions: algorithm runs into problems when probabilities equal zero. To avoid this problem, a small constant (0.0001) is added to each probability in each iteration of the MCMC algorithm, ONLY done for the subject specific probabilities used in the forward algorithm, Multinomial intercept values and group level parameters are left untouched. 
* intercept values for gamma returned for the first iteration of the MCMC algorithm were incorrect: gamma_int_bar[1,], should be transposed before using as.vector. Did not matter for calculations as it is only returned as output to help track the intercept values for gamma over the iteration of the MCMC algorithm. Fixed. 
* fixed all existing problems at https://cran.rstudio.com//web/checks/check_results_mHMMbayes.html -> Escaped LaTeX specials: \& 


# mHMMbayes 0.1.1
Patch release to solve noLD issues (tests without long double on x86_64 Linux system) uncovered by CRAN Package Check Results.

# mHMMbayes 0.1.0 
First (official) version of the package! 
