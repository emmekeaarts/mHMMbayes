
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

You can install the stable mHMMbayes package version directly from CRAN.
Alternatively, you can install the GitHub version that includes the most
recent package developments with:

``` r
# install.packages("devtools")
devtools::install_github("emmekeaarts/mHMMbayes")
```

## Usage

This is a basic example which shows you how to run the model using
categorical example data included with the package. Below we include an
example on continuous data and how to simulate data. For a more
elaborate introduction, see the vignette “tutorial-mhmm” accompanying
the package.

``` r
library(mHMMbayes)

##### Simple 2 state model
# specifying general model properties
m <- 2
n_dep <- 4
q_emiss <- c(3, 2, 3, 2)

# specifying starting values
start_TM <- diag(.8, m)
start_TM[lower.tri(start_TM) | upper.tri(start_TM)] <- .2
start_EM <- list(matrix(c(0.05, 0.90, 0.05,
                          0.90, 0.05, 0.05), byrow = TRUE,
                         nrow = m, ncol = q_emiss[1]), # vocalizing patient
                  matrix(c(0.1, 0.9,
                           0.1, 0.9), byrow = TRUE, nrow = m,
                         ncol = q_emiss[2]), # looking patient
                  matrix(c(0.90, 0.05, 0.05,
                           0.05, 0.90, 0.05), byrow = TRUE,
                         nrow = m, ncol = q_emiss[3]), # vocalizing therapist
                  matrix(c(0.1, 0.9,
                           0.1, 0.9), byrow = TRUE, nrow = m,
                         ncol = q_emiss[4])) # looking therapist

 # Run a model without covariates. 
 # Note that normally, a much higher number of iterations J would be used
 set.seed(23245)
 out_2st <- mHMM(s_data = nonverbal, 
              gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss), 
              start_val = c(list(start_TM), start_EM),
              mcmc = list(J = 11, burn_in = 5))
#> Progress of the Bayesian mHMM algorithm: 
#>   |                                                                              |                                                                      |   0%  |                                                                              |========                                                              |  11%  |                                                                              |================                                                      |  22%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================                                       |  44%  |                                                                              |=======================================                               |  56%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================                |  78%  |                                                                              |==============================================================        |  89%  |                                                                              |======================================================================| 100%
#> Total time elapsed (hh:mm:ss): 00:00:00
 
out_2st
#> Number of subjects: 10 
#> 
#> 11 iterations used in the MCMC algorithm with a burn in of 5 
#> Average Log likelihood over all subjects: -1639.917 
#> Average AIC over all subjects:  3307.834 
#> Average AICc over all subjects: 3308.309 
#> 
#> Number of states used: 2 
#> 
#> Number of dependent variables used: 4 
#> 
#> Type of dependent variable(s): categorical
summary(out_2st)
#> State transition probability matrix 
#>  (at the group level): 
#>  
#>              To state 1 To state 2
#> From state 1      0.915      0.085
#> From state 2      0.075      0.925
#> 
#>  
#> Emission distribution ( categorical ) for each of the dependent variables 
#>  (at the group level): 
#>  
#> $p_vocalizing
#>         Category 1 Category 2 Category 3
#> State 1      0.016      0.970      0.015
#> State 2      0.722      0.124      0.154
#> 
#> $p_looking
#>         Category 1 Category 2
#> State 1      0.220      0.780
#> State 2      0.053      0.947
#> 
#> $t_vocalizing
#>         Category 1 Category 2 Category 3
#> State 1      0.797      0.094      0.110
#> State 2      0.034      0.951      0.015
#> 
#> $t_looking
#>         Category 1 Category 2
#> State 1      0.043      0.957
#> State 2      0.298      0.702

# Run a model including a covariate 
# Here, the covariate (standardized CDI change) predicts the emission 
# distribution for each of the 4 dependent variables:
n_subj <- 10
xx <- rep(list(matrix(1, ncol = 1, nrow = n_subj)), (n_dep + 1))
for(i in 2:(n_dep + 1)){
 xx[[i]] <- cbind(xx[[i]], nonverbal_cov$std_CDI_change)
}
out_2st_c <- mHMM(s_data = nonverbal, xx = xx, 
                 gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss), 
                 start_val = c(list(start_TM), start_EM),
                 mcmc = list(J = 11, burn_in = 5))
#> Progress of the Bayesian mHMM algorithm: 
#>   |                                                                              |                                                                      |   0%  |                                                                              |========                                                              |  11%  |                                                                              |================                                                      |  22%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================                                       |  44%  |                                                                              |=======================================                               |  56%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================                |  78%  |                                                                              |==============================================================        |  89%  |                                                                              |======================================================================| 100%
#> Total time elapsed (hh:mm:ss): 00:00:00
```

This is a basic example which shows you how to run the model on data
that includes continuous (i.e., normally distributed) observations, and
how to simulate (this type of) data.

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
data_cont <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous', gen = list(m = m, n_dep = n_dep),
                  gamma = gamma, emiss_distr = emiss_distr, var_gamma = .1, var_emiss = c(5^2, 0.2^2))

# Specify hyper-prior for the continuous emission distribution
manual_prior_emiss <- prior_emiss_cont(
                        gen = list(m = m, n_dep = n_dep),
                        emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
                                         matrix(c(7, 8, 18), nrow = 1)),
                        emiss_K0 = list(1, 1),
                        emiss_V =  list(rep(5^2, m), rep(0.5^2, m)),
                        emiss_nu = list(1, 1),
                        emiss_a0 = list(rep(1.5, m), rep(1, m)),
                        emiss_b0 = list(rep(20, m), rep(4, m)))

# Run the model on the simulated data:
set.seed(9834)
out_3st_cont_sim <- mHMM(s_data = data_cont$obs,
                         data_distr = 'continuous',
                         gen = list(m = m, n_dep = n_dep),
                         start_val = c(list(gamma), emiss_distr),
                         emiss_hyp_prior = manual_prior_emiss,
                          mcmc = list(J = 1000, burn_in = 200))
#> Progress of the Bayesian mHMM algorithm: 
#>   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%
#> Total time elapsed (hh:mm:ss): 00:00:16

out_3st_cont_sim
#> Number of subjects: 10 
#> 
#> 1000 iterations used in the MCMC algorithm with a burn in of 200 
#> Average Log likelihood over all subjects: -701.7614 
#> Average AIC over all subjects:  1439.523 
#> Average AICc over all subjects: 1447.967 
#> 
#> Number of states used: 3 
#> 
#> Number of dependent variables used: 2 
#> 
#> Type of dependent variable(s): continuous
summary(out_3st_cont_sim)
#> State transition probability matrix 
#>  (at the group level): 
#>  
#>              To state 1 To state 2 To state 3
#> From state 1      0.789      0.081      0.130
#> From state 2      0.222      0.636      0.143
#> From state 3      0.181      0.246      0.574
#> 
#>  
#> Emission distribution ( continuous ) for each of the dependent variables 
#>  (at the group level): 
#>  
#> $`observation 1`
#>            Mean     SD
#> State 1  50.375 11.914
#> State 2  95.105 11.900
#> State 3 154.916 10.212
#> 
#> $`observation 2`
#>           Mean    SD
#> State 1  5.213 2.065
#> State 2  8.968 4.667
#> State 3 19.779 3.162

# obtaining the transition probability matrix gamma and the emission distribution 
# at the group and subject level 
obtain_gamma(out_3st_cont_sim, level = 'group')
#>              To state 1 To state 2 To state 3
#> From state 1      0.789      0.081      0.130
#> From state 2      0.222      0.636      0.143
#> From state 3      0.181      0.246      0.574
obtain_emiss(out_3st_cont_sim, level = 'group')
#> $`observation 1`
#>            Mean     SD
#> State 1  50.375 11.914
#> State 2  95.105 11.900
#> State 3 154.916 10.212
#> 
#> $`observation 2`
#>           Mean    SD
#> State 1  5.213 2.065
#> State 2  8.968 4.667
#> State 3 19.779 3.162

obtain_gamma(out_3st_cont_sim, level = 'subject')
#> $`Subject 1`
#>              To state 1 To state 2 To state 3
#> From state 1      0.800      0.075      0.126
#> From state 2      0.187      0.738      0.076
#> From state 3      0.071      0.158      0.771
#> 
#> $`Subject 2`
#>              To state 1 To state 2 To state 3
#> From state 1      0.855      0.060      0.086
#> From state 2      0.207      0.702      0.091
#> From state 3      0.176      0.312      0.513
#> 
#> $`Subject 3`
#>              To state 1 To state 2 To state 3
#> From state 1      0.671      0.120      0.208
#> From state 2      0.200      0.652      0.149
#> From state 3      0.268      0.184      0.548
#> 
#> $`Subject 4`
#>              To state 1 To state 2 To state 3
#> From state 1      0.819      0.076      0.105
#> From state 2      0.208      0.659      0.133
#> From state 3      0.180      0.267      0.553
#> 
#> $`Subject 5`
#>              To state 1 To state 2 To state 3
#> From state 1      0.869      0.035      0.096
#> From state 2      0.382      0.391      0.228
#> From state 3      0.312      0.161      0.527
#> 
#> $`Subject 6`
#>              To state 1 To state 2 To state 3
#> From state 1      0.856      0.035      0.109
#> From state 2      0.184      0.701      0.116
#> From state 3      0.130      0.446      0.424
#> 
#> $`Subject 7`
#>              To state 1 To state 2 To state 3
#> From state 1      0.797      0.060      0.143
#> From state 2      0.172      0.638      0.190
#> From state 3      0.113      0.228      0.659
#> 
#> $`Subject 8`
#>              To state 1 To state 2 To state 3
#> From state 1      0.835      0.069      0.096
#> From state 2      0.157      0.750      0.093
#> From state 3      0.133      0.138      0.729
#> 
#> $`Subject 9`
#>              To state 1 To state 2 To state 3
#> From state 1      0.732      0.113      0.155
#> From state 2      0.167      0.681      0.152
#> From state 3      0.160      0.310      0.530
#> 
#> $`Subject 10`
#>              To state 1 To state 2 To state 3
#> From state 1      0.860      0.084      0.056
#> From state 2      0.260      0.623      0.117
#> From state 3      0.183      0.272      0.545
obtain_emiss(out_3st_cont_sim, level = 'subject')
#> $`observation 1`
#> $`observation 1`$`Subject 1`
#>            Mean     SD
#> State 1  49.252 11.914
#> State 2  98.191 11.900
#> State 3 159.001 10.212
#> 
#> $`observation 1`$`Subject 2`
#>            Mean     SD
#> State 1  52.368 11.914
#> State 2 109.374 11.900
#> State 3 150.214 10.212
#> 
#> $`observation 1`$`Subject 3`
#>            Mean     SD
#> State 1  54.818 11.914
#> State 2  95.409 11.900
#> State 3 154.159 10.212
#> 
#> $`observation 1`$`Subject 4`
#>            Mean     SD
#> State 1  52.898 11.914
#> State 2  91.227 11.900
#> State 3 153.357 10.212
#> 
#> $`observation 1`$`Subject 5`
#>            Mean     SD
#> State 1  53.216 11.914
#> State 2 101.155 11.900
#> State 3 154.301 10.212
#> 
#> $`observation 1`$`Subject 6`
#>            Mean     SD
#> State 1  51.501 11.914
#> State 2  91.208 11.900
#> State 3 150.539 10.212
#> 
#> $`observation 1`$`Subject 7`
#>            Mean     SD
#> State 1  53.160 11.914
#> State 2  97.857 11.900
#> State 3 154.078 10.212
#> 
#> $`observation 1`$`Subject 8`
#>            Mean     SD
#> State 1  50.592 11.914
#> State 2  97.077 11.900
#> State 3 139.972 10.212
#> 
#> $`observation 1`$`Subject 9`
#>            Mean     SD
#> State 1  55.668 11.914
#> State 2  98.197 11.900
#> State 3 163.877 10.212
#> 
#> $`observation 1`$`Subject 10`
#>            Mean     SD
#> State 1  54.117 11.914
#> State 2  98.876 11.900
#> State 3 152.388 10.212
#> 
#> 
#> $`observation 2`
#> $`observation 2`$`Subject 1`
#>           Mean    SD
#> State 1  4.535 2.065
#> State 2  9.146 4.667
#> State 3 19.726 3.162
#> 
#> $`observation 2`$`Subject 2`
#>           Mean    SD
#> State 1  6.257 2.065
#> State 2  9.061 4.667
#> State 3 20.123 3.162
#> 
#> $`observation 2`$`Subject 3`
#>           Mean    SD
#> State 1  5.030 2.065
#> State 2  9.434 4.667
#> State 3 20.083 3.162
#> 
#> $`observation 2`$`Subject 4`
#>           Mean    SD
#> State 1  4.755 2.065
#> State 2  8.671 4.667
#> State 3 20.273 3.162
#> 
#> $`observation 2`$`Subject 5`
#>           Mean    SD
#> State 1  5.349 2.065
#> State 2  9.162 4.667
#> State 3 19.829 3.162
#> 
#> $`observation 2`$`Subject 6`
#>           Mean    SD
#> State 1  4.843 2.065
#> State 2  9.206 4.667
#> State 3 20.271 3.162
#> 
#> $`observation 2`$`Subject 7`
#>           Mean    SD
#> State 1  4.471 2.065
#> State 2  9.021 4.667
#> State 3 19.458 3.162
#> 
#> $`observation 2`$`Subject 8`
#>           Mean    SD
#> State 1  4.325 2.065
#> State 2  8.164 4.667
#> State 3 19.893 3.162
#> 
#> $`observation 2`$`Subject 9`
#>           Mean    SD
#> State 1  5.455 2.065
#> State 2  9.156 4.667
#> State 3 19.586 3.162
#> 
#> $`observation 2`$`Subject 10`
#>           Mean    SD
#> State 1  5.251 2.065
#> State 2  9.179 4.667
#> State 3 20.127 3.162

# Inferring the most likely state at each point in time
inferred_states <- vit_mHMM(out_3st_cont_sim, data_cont$obs)
#> Please note that the output format is changed from wide to long format to facilitate aditionally returning the state probabilities, see the section 'Value' in the help file for more information.
head(inferred_states)
#>   subj state
#> 1    1     1
#> 2    1     3
#> 3    1     2
#> 4    1     2
#> 5    1     2
#> 6    1     2
```
