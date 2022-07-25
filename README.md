
<!-- README.md is generated from README.Rmd. Please edit that file -->

*Note: this branch is dedicated to development of **using emission
distributions of varying types within the same mHMM model**. When the
development is complete (e.g., documentation is complete, functions
thoroughly tested), the developments will be incorporated in the master
branch, and eventually in the stable R CRAN version of the mHMM
package.*

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

You can install mHMMbayes that incorporates the possibility of using
emission distributions of varying types within the same mHMM model from
github with:

``` r
# install.packages("devtools")
devtools::install_github("emmekeaarts/mHMMbayes@vary-emiss")
```

## Usage

This is a basic example which shows you how to run the model on data
that includes categorical and continuous (i.e., normally distributed)
observations, and how to simulate this type of data.

``` r
library(mHMMbayes)

###### Example on simulated data
# simulating data with both categorical and continous dependent variables

# number of observations per subject, subjects, states and dependent variables
n_t     <- 100
n       <- 10
m       <- 3
n_dep   <- 3

# distribution of each of the dependent variables
data_distr = c("categorical", "continuous", "continuous")

# number of categories for the categorical depenent variable
q_emiss <- c(4, 0, 0)

# transition probabilities between the states
gamma <- matrix(c(0.6, 0.3, 0.1,
                  0.2, 0.7, 0.1,
                  0.1, 0.4, 0.5), nrow = m, byrow = TRUE)

# emission distributions
                    # observation probabilities of each of the 4 categories of the first (categorical) dependent variable over the 3 states
emiss_distr <- list(matrix(c(0.80, 0.13, 0.05, 0.02,
                    0.08, 0.85, 0.05, 0.02,
                    0.05, 0.20, 0.40, 0.35), nrow = m, byrow = TRUE),
                    # mean and variance of second (continuous) dependend in each of the 3 states
                    matrix(c(6.0, 0.5,
                             5.0, 0.5,
                             3.0, 0.5), nrow = m, byrow = TRUE),
                    # mean and variance of third (continous) in each of the 3 states
                    matrix(c(4.0, 1.0,
                           5.0, 1.0,
                           6.5, 1.0), nrow = m, byrow = TRUE))

# amount of variance in the transition probabilities between subjects
var_gamma <- 0.01

# amount of variance in state dependent distributions of the three dependent variables between subjects
var_emiss <- c(0.01, 0.10, 0.10)

# simulate the data
set.seed(3251)
data_vary <- sim_mHMM(n_t = n_t, n = n, data_distr = data_distr, m = m , n_dep = n_dep, q_emiss = q_emiss,
                  gamma = gamma, emiss_distr = emiss_distr, var_gamma = var_gamma, var_emiss = var_emiss)

# Specify hyper-prior for the continuous emission distribution
hyp_pr = emiss_cont_hyp_prior = list(
           emiss_mu0 = list(matrix(c(6,5,3), nrow = 1),
                            matrix(c(4,5,6), nrow = 1)),
           emiss_K0 = list(1,1),
           emiss_nu = list(1,1),
           emiss_V = list(rep(2, 3), rep(2, 3)),
           emiss_a0 = list(rep(1, 3), rep(1, 3)),
           emiss_b0 = list(rep(1, 3), rep(1, 3))
           )

# Run the model on the simulated data:
set.seed(5420)
out_3st_vary_dep <- mHMM_vary(s_data = data_vary$obs,
                    gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                    data_distr = data_distr,
                    start_val = c(list(gamma), emiss_distr),
                    emiss_cont_hyp_prior = hyp_pr,
                    mcmc = list(J = 1000, burn_in = 200))
#> Progress of the Bayesian mHMM algorithm: 
#>   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%
#> Total time elapsed (hh:mm:ss): 00:00:36

out_3st_vary_dep
#> Number of subjects: 10 
#> 
#> 1000 iterations used in the MCMC algorithm with a burn in of 200 
#> Average Log likelihood over all subjects: -403.2717 
#> Average AIC over all subjects: 854.5435 
#> 
#> Number of states used: 3 
#> 
#> Number of dependent variables used: 3 
#> 
#> Distribution(s) of the 3 dependent variable(s) used: 
#>  categorical continuous continuous
summary(out_3st_vary_dep)
#> State transition probability matrix 
#>  (at the group level): 
#>  
#>              To state 1 To state 2 To state 3
#> From state 1      0.540      0.307      0.151
#> From state 2      0.216      0.636      0.143
#> From state 3      0.182      0.349      0.464
#> 
#>  
#> Emission distribution for each of the dependent variables 
#>  (at the group level): 
#>  
#> $`observation 1`
#>         Category 1 Category 2 Category 3 Category 4
#> State 1      0.795      0.152      0.037      0.013
#> State 2      0.056      0.902      0.025      0.011
#> State 3      0.066      0.189      0.418      0.309
#> 
#> $`observation 2`
#>          Mean Variance
#> State 1 5.988    0.774
#> State 2 4.905    0.658
#> State 3 3.085    0.759
#> 
#> $`observation 3`
#>          Mean Variance
#> State 1 4.142    1.012
#> State 2 5.033    1.267
#> State 3 6.185    1.055

# Inferring the most likely state at each point in time
inferred_states <- vit_mHMM_vary(out_3st_vary_dep, data_vary$obs)
head(inferred_states)
#>      Subj_1 Subj_2 Subj_3 Subj_4 Subj_5 Subj_6 Subj_7 Subj_8 Subj_9 Subj_10
#> [1,]      1      1      1      2      3      2      2      2      1       1
#> [2,]      3      1      3      2      2      3      1      1      2       3
#> [3,]      1      3      2      2      3      3      2      1      3       1
#> [4,]      2      3      1      2      2      2      2      1      2       1
#> [5,]      2      1      1      2      2      2      1      1      2       2
#> [6,]      2      1      1      2      3      2      1      1      2       2
```
