
<!-- README.md is generated from README.Rmd. Please edit that file -->

*Note: this branch is dedicated to development of **using continuous
observations (i.e., normally distributed)** within the mHMM model. This
development is near complete. When the development is complete (e.g.,
documentation is complete, functions thoroughly tested), the
developments will be incorporated in the master branch, and eventually
in the stable R CRAN version of the mHMM package.*

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

set.seed(2327)
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
set.seed(9834)
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
#> Average Log likelihood over all subjects: -276.8873 
#> Average AIC over all subjects: 589.7747 
#> 
#> Number of states used: 3 
#> 
#> Number of dependent variables used: 2
summary(out_3st_cont_sim)
#> State transition probability matrix 
#>  (at the group level): 
#>  
#>              To state 1 To state 2 To state 3
#> From state 1      0.755      0.092      0.149
#> From state 2      0.232      0.608      0.156
#> From state 3      0.199      0.259      0.536
#> 
#>  
#> Emission distribution for each of the dependent variables 
#>  (at the group level): 
#>  
#> $`observation 1`
#>           Mean Variance
#> State 1  5.164    1.580
#> State 2  9.454    1.529
#> State 3 15.552    1.160
#> 
#> $`observation 2`
#>          Mean Variance
#> State 1 0.516    0.133
#> State 2 0.890    0.234
#> State 3 2.019    0.152

# obtaining the transition probability matrix gamma and the emission ditribution 
# at the group and subject level 
obtain_gamma(out_3st_cont_sim, level = 'group')
#>              To state 1 To state 2 To state 3
#> From state 1      0.755      0.092      0.149
#> From state 2      0.232      0.608      0.156
#> From state 3      0.199      0.259      0.536
obtain_emiss(out_3st_cont_sim, level = 'subject')
#> $`observation 1`
#> $`observation 1`$`Subject 1`
#>           Mean Variance
#> State 1  4.843    1.580
#> State 2  9.798    1.529
#> State 3 16.234    1.160
#> 
#> $`observation 1`$`Subject 2`
#>           Mean Variance
#> State 1  5.298    1.580
#> State 2 11.209    1.529
#> State 3 14.965    1.160
#> 
#> $`observation 1`$`Subject 3`
#>           Mean Variance
#> State 1  5.586    1.580
#> State 2  9.603    1.529
#> State 3 15.503    1.160
#> 
#> $`observation 1`$`Subject 4`
#>           Mean Variance
#> State 1  5.425    1.580
#> State 2  8.746    1.529
#> State 3 15.436    1.160
#> 
#> $`observation 1`$`Subject 5`
#>           Mean Variance
#> State 1  5.521    1.580
#> State 2 10.122    1.529
#> State 3 15.573    1.160
#> 
#> $`observation 1`$`Subject 6`
#>           Mean Variance
#> State 1  5.222    1.580
#> State 2  8.762    1.529
#> State 3 15.003    1.160
#> 
#> $`observation 1`$`Subject 7`
#>           Mean Variance
#> State 1  5.458    1.580
#> State 2  9.716    1.529
#> State 3 15.470    1.160
#> 
#> $`observation 1`$`Subject 8`
#>           Mean Variance
#> State 1  5.127    1.580
#> State 2  9.687    1.529
#> State 3 13.646    1.160
#> 
#> $`observation 1`$`Subject 9`
#>           Mean Variance
#> State 1  5.822    1.580
#> State 2  9.662    1.529
#> State 3 16.903    1.160
#> 
#> $`observation 1`$`Subject 10`
#>           Mean Variance
#> State 1  5.621    1.580
#> State 2  9.954    1.529
#> State 3 15.327    1.160
#> 
#> 
#> $`observation 2`
#> $`observation 2`$`Subject 1`
#>          Mean Variance
#> State 1 0.318    0.133
#> State 2 0.871    0.234
#> State 3 2.030    0.152
#> 
#> $`observation 2`$`Subject 2`
#>          Mean Variance
#> State 1 0.872    0.133
#> State 2 0.806    0.234
#> State 3 2.038    0.152
#> 
#> $`observation 2`$`Subject 3`
#>          Mean Variance
#> State 1 0.474    0.133
#> State 2 1.014    0.234
#> State 3 1.942    0.152
#> 
#> $`observation 2`$`Subject 4`
#>          Mean Variance
#> State 1 0.424    0.133
#> State 2 0.854    0.234
#> State 3 2.113    0.152
#> 
#> $`observation 2`$`Subject 5`
#>          Mean Variance
#> State 1 0.615    0.133
#> State 2 0.952    0.234
#> State 3 2.006    0.152
#> 
#> $`observation 2`$`Subject 6`
#>          Mean Variance
#> State 1 0.445    0.133
#> State 2 0.937    0.234
#> State 3 2.119    0.152
#> 
#> $`observation 2`$`Subject 7`
#>          Mean Variance
#> State 1 0.384    0.133
#> State 2 0.957    0.234
#> State 3 2.035    0.152
#> 
#> $`observation 2`$`Subject 8`
#>          Mean Variance
#> State 1 0.393    0.133
#> State 2 0.721    0.234
#> State 3 2.107    0.152
#> 
#> $`observation 2`$`Subject 9`
#>          Mean Variance
#> State 1 0.545    0.133
#> State 2 0.861    0.234
#> State 3 1.923    0.152
#> 
#> $`observation 2`$`Subject 10`
#>          Mean Variance
#> State 1 0.590    0.133
#> State 2 0.964    0.234
#> State 3 2.076    0.152

# Inferring the most likely state at each point in time
inferred_states <- vit_mHMM_cont(out_3st_cont_sim, data_cont$obs)
head(inferred_states)
#>      Subj_1 Subj_2 Subj_3 Subj_4 Subj_5 Subj_6 Subj_7 Subj_8 Subj_9 Subj_10
#> [1,]      1      1      1      2      1      1      3      2      1       1
#> [2,]      3      1      1      1      1      1      3      2      1       1
#> [3,]      2      1      1      1      1      1      3      3      1       1
#> [4,]      2      2      1      1      1      1      3      3      1       1
#> [5,]      2      2      3      1      1      1      3      3      1       1
#> [6,]      2      1      3      1      1      1      3      3      1       2
```
