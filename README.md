
<!-- README.md is generated from README.Rmd. Please edit that file -->

\*Note: this branch is dedicated to development of \_\_using continuous
observations (i.e., normally distributed) within the mHMM model. This
development is near complete. When the development is complete (e.g.,
documentation is complete, functions thoroughly tested), the
developments will be incorporated in the master branch, and eventually
in the stable R CRAN version of the mHMM package.\*

# mHMMbayes

With the package mHMMbayes you can fit multilevel hidden Markov models.
The multilevel hidden Markov model (HMM) is a generalization of the
well-known hidden Markov model, tailored to accommodate (intense)
longitudinal data of multiple individuals simultaneously. Using a
multilevel framework, we allow for heterogeneity in the model parameters
(transition probability matrix and conditional distribution), while
estimating one overall HMM. The model has a great potential of
application in many fields, such as the social sciences and medicine.
The model can be fitted on multivariate data with categorical and/or
continuous (i.e., normally distributed) observations, and include
individual level covariates (allowing for e.g., group comparisons on
model parameters). Parameters are estimated using Bayesian estimation
utilizing the forward-backward recursion within a hybrid Metropolis
within Gibbs sampler. The package also includes various options for
model visualization, a function to simulate data and a function to
obtain the most likely hidden state sequence for each individual using
the Viterbi algorithm.

Please do not hesitate to contact me if you have any questions regarding
the package.

## Installation

You can install mHMMbayes that incorporates the possibility of modelling
data with continuous observations from github with:

``` r
# install.packages("devtools")
devtools::install_github("emmekeaarts/mHMMbayes@continuous-emiss")
```

## Usage

This is a basic example which shows you how to run the model on data
that includes continuous (i.e., normally distributed) observations, and
how to simulate this type of data.

``` r
library(mHMMbayes)

###### Example on simulated data
# simulating multivariate continuous data
n_t     <- 100
n       <- 10
m       <- 3
n_dep   <- 2

gamma   <- matrix(c(0.8, 0.1, 0.1,
                    0.2, 0.7, 0.1,
                    0.2, 0.2, 0.6), ncol = m, byrow = TRUE)

emiss_distr <- list(matrix(c( 5, 1,
                              10, 1,
                              15, 1), nrow = m, byrow = TRUE),
                    matrix(c(0.5, 0.1,
                             1.0, 0.2,
                             2.0, 0.1), nrow = m, byrow = TRUE))

data_cont <- sim_mHMM(n_t = n_t, n = n, m = m, n_dep = n_dep, data_distr = 'continuous',
                  gamma = gamma, emiss_distr = emiss_distr, var_gamma = .1, var_emiss = c(.5, 0.01))

# Specify hyper-prior for the continuous emission distribution
hyp_pr <- list(
               emiss_mu0 = list(matrix(c(3,7,17), nrow = 1), matrix(c(0.7, 0.8, 1.8), nrow = 1)),
               emiss_K0  = list(1, 1),
               emiss_nu  = list(1, 1),
               emiss_V   = list(rep(2, m), rep(1, m)),
               emiss_a0  = list(rep(1, m), rep(1, m)),
               emiss_b0  = list(rep(1, m), rep(1, m)))

# Run the model on the simulated data:
out_3st_cont_sim <- mHMM_cont(s_data = data_cont$obs,
                    gen = list(m = m, n_dep = n_dep),
                    start_val = c(list(gamma), emiss_distr),
                    emiss_hyp_prior = hyp_pr,
                    mcmc = list(J = 1000, burn_in = 200))
#> Progress of the Bayesian mHMM algorithm: 
#>   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%
#> Total time elapsed (hh:mm:ss): 00:00:25

out_3st_cont_sim
#> Number of subjects: 10 
#> 
#> 1000 iterations used in the MCMC algorithm with a burn in of 200 
#> Average Log likelihood over all subjects: -286.2929 
#> Average AIC over all subjects: 608.5859 
#> 
#> Number of states used: 3 
#> 
#> Number of dependent variables used: 2
summary(out_3st_cont_sim)
#> State transition probability matrix 
#>  (at the group level): 
#>  
#>              To state 1 To state 2 To state 3
#> From state 1      0.729      0.124      0.143
#> From state 2      0.225      0.601      0.171
#> From state 3      0.239      0.227      0.532
#> 
#>  
#> Emission distribution for each of the dependent variables 
#>  (at the group level): 
#>  
#> $`observation 1`
#>           Mean Variance
#> State 1  4.842    1.314
#> State 2  9.402    1.447
#> State 3 15.435    1.555
#> 
#> $`observation 2`
#>          Mean Variance
#> State 1 0.448    0.143
#> State 2 0.932    0.255
#> State 3 1.992    0.156

# obtaining the transition probability matrix gamma and the emission ditribution 
# at the group and subject level 
obtain_gamma(out_3st_cont_sim, level = 'group')
#>              To state 1 To state 2 To state 3
#> From state 1      0.729      0.124      0.143
#> From state 2      0.225      0.601      0.171
#> From state 3      0.239      0.227      0.532
obtain_emiss(out_3st_cont_sim, level = 'subject')
#> $`observation 1`
#> $`observation 1`$`Subject 1`
#>           Mean Variance
#> State 1  5.275    1.314
#> State 2 10.497    1.447
#> State 3 15.211    1.555
#> 
#> $`observation 1`$`Subject 2`
#>           Mean Variance
#> State 1  5.440    1.314
#> State 2  8.936    1.447
#> State 3 16.373    1.555
#> 
#> $`observation 1`$`Subject 3`
#>           Mean Variance
#> State 1  4.649    1.314
#> State 2  9.871    1.447
#> State 3 15.153    1.555
#> 
#> $`observation 1`$`Subject 4`
#>           Mean Variance
#> State 1  5.092    1.314
#> State 2  9.380    1.447
#> State 3 15.714    1.555
#> 
#> $`observation 1`$`Subject 5`
#>           Mean Variance
#> State 1  5.800    1.314
#> State 2  9.480    1.447
#> State 3 15.259    1.555
#> 
#> $`observation 1`$`Subject 6`
#>           Mean Variance
#> State 1  4.912    1.314
#> State 2  9.242    1.447
#> State 3 14.182    1.555
#> 
#> $`observation 1`$`Subject 7`
#>           Mean Variance
#> State 1  5.282    1.314
#> State 2  9.382    1.447
#> State 3 14.506    1.555
#> 
#> $`observation 1`$`Subject 8`
#>           Mean Variance
#> State 1  5.324    1.314
#> State 2  9.694    1.447
#> State 3 17.132    1.555
#> 
#> $`observation 1`$`Subject 9`
#>           Mean Variance
#> State 1  3.871    1.314
#> State 2 10.121    1.447
#> State 3 14.221    1.555
#> 
#> $`observation 1`$`Subject 10`
#>           Mean Variance
#> State 1  4.656    1.314
#> State 2  9.760    1.447
#> State 3 14.916    1.555
#> 
#> 
#> $`observation 2`
#> $`observation 2`$`Subject 1`
#>          Mean Variance
#> State 1 0.507    0.143
#> State 2 0.922    0.255
#> State 3 1.827    0.156
#> 
#> $`observation 2`$`Subject 2`
#>          Mean Variance
#> State 1 0.542    0.143
#> State 2 0.770    0.255
#> State 3 2.216    0.156
#> 
#> $`observation 2`$`Subject 3`
#>          Mean Variance
#> State 1 0.327    0.143
#> State 2 0.813    0.255
#> State 3 1.935    0.156
#> 
#> $`observation 2`$`Subject 4`
#>          Mean Variance
#> State 1 0.354    0.143
#> State 2 0.954    0.255
#> State 3 1.873    0.156
#> 
#> $`observation 2`$`Subject 5`
#>          Mean Variance
#> State 1 0.347    0.143
#> State 2 0.891    0.255
#> State 3 1.984    0.156
#> 
#> $`observation 2`$`Subject 6`
#>          Mean Variance
#> State 1 0.431    0.143
#> State 2 0.933    0.255
#> State 3 1.869    0.156
#> 
#> $`observation 2`$`Subject 7`
#>          Mean Variance
#> State 1 0.430    0.143
#> State 2 1.130    0.255
#> State 3 2.044    0.156
#> 
#> $`observation 2`$`Subject 8`
#>          Mean Variance
#> State 1 0.333    0.143
#> State 2 1.008    0.255
#> State 3 2.128    0.156
#> 
#> $`observation 2`$`Subject 9`
#>          Mean Variance
#> State 1 0.674    0.143
#> State 2 0.933    0.255
#> State 3 2.098    0.156
#> 
#> $`observation 2`$`Subject 10`
#>          Mean Variance
#> State 1 0.257    0.143
#> State 2 0.987    0.255
#> State 3 2.124    0.156

# Inferring the most likely state at each point in time
inferred_states <- vit_mHMM_cont(out_3st_cont_sim, data_cont$obs)
head(inferred_states)
#>      Subj_1 Subj_2 Subj_3 Subj_4 Subj_5 Subj_6 Subj_7 Subj_8 Subj_9 Subj_10
#> [1,]      2      1      1      1      3      3      1      1      1       3
#> [2,]      2      1      1      1      3      3      1      1      2       3
#> [3,]      3      2      1      1      1      3      2      1      2       3
#> [4,]      2      2      2      1      1      1      3      1      2       3
#> [5,]      3      2      2      2      1      1      3      2      2       1
#> [6,]      3      2      2      2      1      1      3      2      2       1
```
