
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/emmekeaarts/mHMMbayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/emmekeaarts/mHMMbayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# mHMMbayes

With the package mHMMbayes you can fit multilevel hidden Markov models.
The multilevel hidden Markov model (HMM) is a generalization of the
well-known hidden Markov model, tailored to accommodate (intense)
longitudinal data of multiple individuals simultaneously. Using a
multilevel framework, we allow for heterogeneity in the model parameters
(transition probability matrix and conditional distribution), while
estimating one overall HMM. The model has a great potential of
application in many fields, such as the social sciences and medicine.
The model can be fitted on multivariate data with categorical,
continuous (i.e., normally distributed), or count (i.e., Poisson
distributed) observations, and include individual level covariates
(allowing for e.g., group comparisons on model parameters). Parameters
are estimated using Bayesian estimation utilizing the forward-backward
recursion within a hybrid Metropolis within Gibbs sampler. Missing data
(NA) in the dependent variables is accommodated assuming MAR. The
package also includes various options for model visualization, a
function to simulate data and a function to obtain the most likely
hidden state sequence for each individual using the Viterbi algorithm.

Please do not hesitate to contact me if you have any questions regarding
the package.

## Installation

You can install mHMMbayes that incorporates the possibility of fitting a
1-state model from github with:

``` r
# install.packages("devtools")
devtools::install_github("emmekeaarts/mHMMbayes@one-state")
```

## Usage

This is a basic example which shows you how to run the model using
simulated continuous data with one state. For a general, more elaborate
introduction, see the vignette “tutorial-mhmm” accompanying the package.

``` r
library(mHMMbayes)

###### Example on simulated data
# simulating multivariate continuous data
n_t     <- 100
n       <- 10
m       <- 1
n_dep   <- 2

gamma   <- matrix(c(1), ncol = m, byrow = TRUE)

emiss_distr <- list(matrix(c( 50, 10), nrow = m, byrow = TRUE),
                    matrix(c(5, 2), nrow = m, byrow = TRUE))

set.seed(2327)
data_cont <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous', gen = list(m = m, n_dep = n_dep),
                  gamma = gamma, emiss_distr = emiss_distr, var_gamma = .1, var_emiss = c(5^2, 0.2^2))
#> Warning in sim_mHMM(n_t = n_t, n = n, data_distr = "continuous", gen = list(m =
#> m, : If the number of states m is set to 1, variance between subjects in the
#> transition probability matrix is ignored, and cannot be explained by
#> covariates.

# Specify hyper-prior for the continuous emission distribution
manual_prior_emiss <- prior_emiss_cont(gen = list(m = m, n_dep = n_dep),
                                       emiss_mu0 = list(matrix(c(30), nrow = 1),
                                                        matrix(c(7), nrow = 1)),
                                       emiss_K0 = list(1, 1),
                                       emiss_V =  list(rep(5^2, m), rep(0.5^2, m)),
                                       emiss_nu = list(1, 1),
                                       emiss_a0 = list(rep(1.5, m), rep(1, m)),
                                       emiss_b0 = list(rep(20, m), rep(4, m)))

# Run the model on the simulated data:
set.seed(9834)
out_1st_cont_sim <- mHMM(s_data = data_cont$obs,
                         data_distr = 'continuous',
                         gen = list(m = m, n_dep = n_dep),
                         start_val = c(list(gamma), emiss_distr),
                         emiss_hyp_prior = manual_prior_emiss,
                          mcmc = list(J = 1000, burn_in = 200))
#> Progress of the Bayesian mHMM algorithm: 
#>   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%
#> Total time elapsed (hh:mm:ss): 00:00:00

out_1st_cont_sim
#> Number of subjects: 10 
#> 
#> 1000 iterations used in the MCMC algorithm with a burn in of 200 
#> Average Log likelihood over all subjects: -598.6083 
#> Average AIC over all subjects:  1205.217 
#> Average AICc over all subjects: 1205.638 
#> 
#> Number of states used: 1 
#> 
#> Number of dependent variables used: 2 
#> 
#> Type of dependent variable(s): continuous
summary(out_1st_cont_sim)
#> State transition probability matrix 
#>  (at the group level): 
#>  
#>              To state 1
#> From state 1          1
#> 
#>  
#> Emission distribution ( continuous ) for each of the dependent variables 
#>  (at the group level): 
#>  
#> $`observation 1`
#>           Mean     SD
#> State 1 45.093 11.489
#> 
#> $`observation 2`
#>          Mean    SD
#> State 1 5.163 2.131

# obtaining the emission distribution 
# at the group and subject level 
obtain_emiss(out_1st_cont_sim, level = 'group')
#> $`observation 1`
#>           Mean     SD
#> State 1 45.093 11.489
#> 
#> $`observation 2`
#>          Mean    SD
#> State 1 5.163 2.131
obtain_emiss(out_1st_cont_sim, level = 'subject')
#> $`observation 1`
#> $`observation 1`$`Subject 1`
#>          Mean     SD
#> State 1 43.79 11.489
#> 
#> $`observation 1`$`Subject 2`
#>           Mean     SD
#> State 1 46.645 11.489
#> 
#> $`observation 1`$`Subject 3`
#>           Mean     SD
#> State 1 46.747 11.489
#> 
#> $`observation 1`$`Subject 4`
#>           Mean     SD
#> State 1 39.139 11.489
#> 
#> $`observation 1`$`Subject 5`
#>           Mean     SD
#> State 1 45.154 11.489
#> 
#> $`observation 1`$`Subject 6`
#>           Mean     SD
#> State 1 59.812 11.489
#> 
#> $`observation 1`$`Subject 7`
#>           Mean     SD
#> State 1 47.716 11.489
#> 
#> $`observation 1`$`Subject 8`
#>           Mean     SD
#> State 1 45.042 11.489
#> 
#> $`observation 1`$`Subject 9`
#>           Mean     SD
#> State 1 47.048 11.489
#> 
#> $`observation 1`$`Subject 10`
#>           Mean     SD
#> State 1 45.065 11.489
#> 
#> 
#> $`observation 2`
#> $`observation 2`$`Subject 1`
#>          Mean    SD
#> State 1 4.969 2.131
#> 
#> $`observation 2`$`Subject 2`
#>          Mean    SD
#> State 1 5.306 2.131
#> 
#> $`observation 2`$`Subject 3`
#>          Mean    SD
#> State 1 5.053 2.131
#> 
#> $`observation 2`$`Subject 4`
#>          Mean    SD
#> State 1 4.608 2.131
#> 
#> $`observation 2`$`Subject 5`
#>          Mean    SD
#> State 1 4.689 2.131
#> 
#> $`observation 2`$`Subject 6`
#>          Mean    SD
#> State 1 4.936 2.131
#> 
#> $`observation 2`$`Subject 7`
#>          Mean    SD
#> State 1 4.661 2.131
#> 
#> $`observation 2`$`Subject 8`
#>          Mean    SD
#> State 1 4.773 2.131
#> 
#> $`observation 2`$`Subject 9`
#>          Mean    SD
#> State 1 5.402 2.131
#> 
#> $`observation 2`$`Subject 10`
#>          Mean    SD
#> State 1 5.283 2.131

# note that as the transition probability matrix always equals 1, 
# as such running the function obtain_gamma() function on output of a 1-state model will not run (and will throw an explanation)

# note that as all states over time equal 1, inferring the most likely state at each point in time is obsolete. 
# as such, running the vit_mHMM() function on output of a 1-state model will not run (and will throw an explanation)
```
