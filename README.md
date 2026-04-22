
<!-- README.md is generated from README.Rmd. Please edit that file -->

*Note: this branch is dedicated to development of **fitting a model with
FIXED GAUSSIAN emission distributions** using the `mHMM_f()` function,
creating an object of class `mHMM_f`. In this branch, the usual function
`mHMM()` is not included. All S3 methods are tailored to the `mHMM_f()`
function.*

# mHMMbayes

With the package mHMMbayes you can fit multilevel hidden Markov models.
The multilevel hidden Markov model (HMM) is a generalization of the
well-known hidden Markov model, tailored to accommodate (intense)
longitudinal data of multiple individuals simultaneously. Using a
multilevel framework, we allow for heterogeneity in the model parameters
(transition probability matrix), while estimating one overall HMM. The
model has a great potential of application in many fields, such as the
social sciences and medicine. The model (in this specific branch) can be
fitted on multivariate data with continuous (i.e., normally distributed)
observations, and include individual level covariates (allowing for
e.g., group comparisons on model parameters). Parameters are estimated
using Bayesian estimation utilizing the forward-backward recursion
within a hybrid Metropolis within Gibbs sampler. Missing data (NA) in
the dependent variables is accommodated assuming MAR. The package also
includes various options for model visualization, a function to simulate
data and a function to obtain the most likely hidden state sequence for
each individual using the Viterbi algorithm.

Please do not hesitate to contact me if you have any questions regarding
the package.

## Installation

You can install mHMMbayes that accommodates fixed over subjects Gaussian
emission distributions from github with:

``` r
# install.packages("devtools")
devtools::install_github("emmekeaarts/mHMMbayes@fixed-emiss")
```

## Usage

This is a basic example which shows you how to run the model using
continuous data and how to simulate data.

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

emiss_distr <- list(matrix(c( 50, 10,
                              100, 10,
                              150, 10), nrow = m, byrow = TRUE),
                    matrix(c(5, 2,
                             10, 5,
                             20, 3), nrow = m, byrow = TRUE))

set.seed(2327)
data_cont <- sim_mHMM_f(n_t = n_t, n = n, data_distr = 'continuous', gen = list(m = m, n_dep = n_dep),
                  gamma = gamma, emiss_distr = emiss_distr, var_gamma = .1)

# Specify hyper-prior for the continuous emission distribution
manual_prior_emiss <- prior_emiss_cont_f(
                        gen = list(m = m, n_dep = n_dep),
                        emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
                                         matrix(c(7, 8, 18), nrow = 1)),
                        emiss_K0 = list(1, 1),
                        emiss_V =  list(rep(5^2, m), rep(0.5^2, m)),
                        emiss_nu = list(1, 1))

# Run the model on the simulated data:
set.seed(9834)
out_3st_cont_sim <- mHMM_f(s_data = data_cont$obs,
                         data_distr = 'continuous',
                         gen = list(m = m, n_dep = n_dep),
                         start_val = c(list(gamma), emiss_distr),
                         emiss_hyp_prior = manual_prior_emiss,
                          mcmc = list(J = 1000, burn_in = 200))
#> Progress of the Bayesian mHMM algorithm: 
#>   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%
#> Total time elapsed (hh:mm:ss): 00:00:15

out_3st_cont_sim
#> Number of subjects: 10 
#> 
#> 1000 iterations used in the MCMC algorithm with a burn in of 200 
#> Average Log likelihood over all subjects: -678.0322 
#> Average AIC over all subjects:  1392.064 
#> Average AICc over all subjects: 1400.509 
#> 
#> Number of states used: 3 
#> 
#> Number of dependent variables used: 2 
#> 
#> Type of dependent variable(s): continuous (fixed emission)
summary(out_3st_cont_sim)
#> State transition probability matrix 
#>  (at the group level): 
#>  
#>              To state 1 To state 2 To state 3
#> From state 1      0.828      0.089      0.084
#> From state 2      0.188      0.706      0.106
#> From state 3      0.202      0.206      0.591
#> 
#>  
#> Emission distribution ( continuous ) for each of the dependent variables 
#>  (at the group level): 
#>  
#> $`observation 1`
#>            Mean     SD
#> State 1  50.483 10.162
#> State 2  99.950  9.789
#> State 3 149.131  9.167
#> 
#> $`observation 2`
#>           Mean    SD
#> State 1  4.989 1.932
#> State 2 10.376 5.140
#> State 3 19.676 3.067

# obtaining the transition probability matrix gamma and the emission distribution 
# at the group level 
obtain_gamma(out_3st_cont_sim, level = 'group')
#>              To state 1 To state 2 To state 3
#> From state 1      0.828      0.089      0.084
#> From state 2      0.188      0.706      0.106
#> From state 3      0.202      0.206      0.591
obtain_emiss_f(out_3st_cont_sim)       # level is not specified here as only group level parameters are available. 
#> $`observation 1`
#>            Mean     SD
#> State 1  50.483 10.162
#> State 2  99.950  9.789
#> State 3 149.131  9.167
#> 
#> $`observation 2`
#>           Mean    SD
#> State 1  4.989 1.932
#> State 2 10.376 5.140
#> State 3 19.676 3.067

# and ONLY for the transition probabilities at the subject level
obtain_gamma(out_3st_cont_sim, level = 'subject')
#> $`Subject 1`
#>              To state 1 To state 2 To state 3
#> From state 1      0.941      0.041      0.018
#> From state 2      0.090      0.874      0.035
#> From state 3      0.166      0.249      0.585
#> 
#> $`Subject 2`
#>              To state 1 To state 2 To state 3
#> From state 1      0.779      0.084      0.137
#> From state 2      0.167      0.794      0.038
#> From state 3      0.188      0.244      0.569
#> 
#> $`Subject 3`
#>              To state 1 To state 2 To state 3
#> From state 1      0.812      0.089      0.100
#> From state 2      0.207      0.589      0.204
#> From state 3      0.205      0.166      0.629
#> 
#> $`Subject 4`
#>              To state 1 To state 2 To state 3
#> From state 1      0.826      0.110      0.064
#> From state 2      0.192      0.620      0.188
#> From state 3      0.331      0.248      0.422
#> 
#> $`Subject 5`
#>              To state 1 To state 2 To state 3
#> From state 1      0.837      0.058      0.105
#> From state 2      0.283      0.533      0.184
#> From state 3      0.104      0.081      0.815
#> 
#> $`Subject 6`
#>              To state 1 To state 2 To state 3
#> From state 1      0.864      0.050      0.086
#> From state 2      0.086      0.850      0.064
#> From state 3      0.191      0.099      0.710
#> 
#> $`Subject 7`
#>              To state 1 To state 2 To state 3
#> From state 1      0.770      0.104      0.126
#> From state 2      0.156      0.779      0.064
#> From state 3      0.114      0.125      0.761
#> 
#> $`Subject 8`
#>              To state 1 To state 2 To state 3
#> From state 1      0.880      0.058      0.062
#> From state 2      0.218      0.655      0.127
#> From state 3      0.142      0.316      0.542
#> 
#> $`Subject 9`
#>              To state 1 To state 2 To state 3
#> From state 1      0.842      0.100      0.058
#> From state 2      0.175      0.778      0.047
#> From state 3      0.265      0.285      0.450
#> 
#> $`Subject 10`
#>              To state 1 To state 2 To state 3
#> From state 1      0.877      0.085      0.039
#> From state 2      0.175      0.699      0.127
#> From state 3      0.233      0.241      0.526

# Inferring the most likely state at each point in time
inferred_states <- vit_mHMM_f(out_3st_cont_sim, data_cont$obs)
#> Please note that the output format is changed from wide to long format to facilitate aditionally returning the state probabilities, see the section 'Value' in the help file for more information.
head(inferred_states)
#>   subj state
#> 1    1     1
#> 2    1     1
#> 3    1     1
#> 4    1     1
#> 5    1     1
#> 6    1     2
```
