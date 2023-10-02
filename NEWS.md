# mHMMbayes 0.3.0

## Accommodating continous data (i.e., Normally distributed data)
A major improvement in this release is the possibility to include continuous data in `mHMM()`. Currently, the user can model data composed of either categorical data OR continuous data (so a mix of different types of emission distributions is not possible within the stable CRAN version). As such, the following changes are implemented:

* `mHMM()` now includes the input parameter `data_distr`, where the user can specify whether the input data contains categorical or continuous data. Defaults to `data_distr = 'categorical'`.
* `sim_mHMM()` allows the simulation of continuous data, which is facilitated again by the input parameter `data_distr`, where the user can specify whether one wants to simulate categorical or continuous data. Defaults to `data_distr = 'categorical'`.
* Functions that utilize `mHMM()` output objects as input such as `obtain_emiss()`, `vit_mHMM()`, and S3 methods as `print()`, `summary()`, and `plot()` automatically detect whether the output object relates to a multilevel HMM fitted to categorical or continuous data, and adjusts it's processing methods accordingly.  

Also, a new function, `prior_emiss_cont()`, is introduced which enables the specification of hyper prior parameters when using continuous input data. 

## Missing data (NA) in the dependent variable(s) 
New is also the accommodation of missing values (`NA`) in the dependent input variable(s). Missingness is assumed Missing at Random (MAR), so that the missingness mechanism is independent of the missing data and the hidden states given the observed data and model parameters. This means that missing observations are assumed equally likely in each of the states. In our approach, hidden state probabilities are inferred for missing observations (thus only based on the transition probability matrix gamma), but missing observations themselves are not directly imputed. 

To accommodate missing values, the forward algorithm implemented in C++ was slightly adjusted, and state dependent observations (on which the parameter estimates of the emission distribution(s) are based) are selected such that missing values are omitted. 

## Returned output by `mHMM()`
The `mHMM()` output component `PD_subj` was modified to facilitate the inclusion of both categorical and continuous input data. Before, `PD_subj` was a list containing one matrix per subject containing all subject level output parameters over the iterations of the MCMC sampler. Now, `PD_subj`is a list containing one list per subject with the elements `trans_prob`, `cat_emiss` or `cont_emiss` in case of categorical or continuous observations, respectively, and `log_likl`, providing the subject parameter estimates over the iterations of the MCMC sampler. `trans_prob` relates to the transition probabilities gamma, `cat_emiss` to the categorical emission distribution (emission probabilities), `cont_emiss` to the continuous emission distributions (subsequently the the emission means and the (fixed over subjects) emission standard deviation), and `log_likl` to the log likelihood over the MCMC iterations. 

A detailed error message is displayed when trying to post-process mHMM objects created with an earlier version of the package. 

## Extra checks input data within `mHMM()`
Several extra checks have been implemented in `mHMM()`. Specifically, checking for: 

* The correct specification of categorical input, where values are allowed to range from 1 to the number of categories observable within a variable.
* The inclusion of zero values in the start values of the transition prabability matrix gamma or the categorical emission probabilities, as this can lead to problems in the forward algorithm.
* Possible reasons for a fatal error in the forward algorithm in the first (starting values that do not sufficiently support the range of observed values) or subsequent (hyper-paramter values that result in hyper distributions that do not sufficiently support the range of observed values) iterations. 
* The correct dimensions for the starting values specified via `starting_val`.

## Other minor (quite technical) improvements:
* The appropriate PKGNAME-package \alias as per "Documenting packages" is now included. 
* The minimum supported version of R is now 3.6 (instead of 3.5), as the package now requires compilation with at least C++11 (and R 3.6, defaults to compiling packages with the C++11 standard). 
* First, to the state transitions of the sampled state sequence in the MCMC sampler, and to state dependent categorical observations, a sequence of 1:m and 1:q_emiss[k] was added to ensure that all possible outcomes were observed at least once to avoid estimation problems (e.g., when a certain state was not sampled at all within an iteration, the sampled state sequence would equal 1:m to avoid 'empty' state transition observations). However, this resulted in bias, and this approach is now completely omitted. Parameter estimates are now only based on the sampled state sequence and sampled state dependent observations, in combination with the prior distribution. This means that if a certain state is not observed at all, parameter estimation solely depends on the prior distribution. 


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
