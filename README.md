
<!-- README.md is generated from README.Rmd. Please edit that file -->

*Note: In contrast to the current CRAN version of mHMMbayes, this GitHub
version allows for **using continuous observations (i.e., normally
distributed)** and missing data (NA) within the mHMM model. When this
new version is even more thoroughly tested, the developments will be
incorporated in the stable R CRAN version of the mHMM package.*

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
within Gibbs sampler. Missing data (NA) in the dependent variables is
accommodated assuming MAR. The package also includes various options for
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
#> Total time elapsed (hh:mm:ss): 00:00:01
 
out_2st
#> Number of subjects: 10 
#> 
#> 11 iterations used in the MCMC algorithm with a burn in of 5 
#> Average Log likelihood over all subjects: -1639.917 
#> Average AIC over all subjects: 3307.834 
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
#> Total time elapsed (hh:mm:ss): 00:00:01
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
                  gamma = gamma, emiss_distr = emiss_distr, var_gamma = .1, var_emiss = c(.5, 0.01))

# Specify hyper-prior for the continuous emission distribution
manual_prior_emiss <- prior_emiss_cont(
                        gen = list(m = m, n_dep = n_dep),
                        emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
                                         matrix(c(7, 8, 18), nrow = 1)),
                        emiss_K0 = list(1, 1),
                        emiss_V =  list(rep(100, m), rep(25, m)),
                        emiss_nu = list(1, 1),
                        emiss_a0 = list(rep(1, m), rep(1, m)),
                        emiss_b0 = list(rep(1, m), rep(1, m)))

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
#> Total time elapsed (hh:mm:ss): 00:00:17

out_3st_cont_sim
#> Number of subjects: 10 
#> 
#> 1000 iterations used in the MCMC algorithm with a burn in of 200 
#> Average Log likelihood over all subjects: -707.9179 
#> Average AIC over all subjects: 1451.836 
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
#> From state 1      0.795      0.076      0.128
#> From state 2      0.218      0.645      0.137
#> From state 3      0.179      0.252      0.570
#> 
#>  
#> Emission distribution ( continuous ) for each of the dependent variables 
#>  (at the group level): 
#>  
#> $`observation 1`
#>            Mean     SD
#> State 1  48.543 11.778
#> State 2  95.767 11.797
#> State 3 153.027 10.156
#> 
#> $`observation 2`
#>           Mean    SD
#> State 1  5.227 2.226
#> State 2  9.218 4.954
#> State 3 19.812 3.389

# obtaining the transition probability matrix gamma and the emission ditribution 
# at the group and subject level 
obtain_gamma(out_3st_cont_sim, level = 'group')
#>              To state 1 To state 2 To state 3
#> From state 1      0.795      0.076      0.128
#> From state 2      0.218      0.645      0.137
#> From state 3      0.179      0.252      0.570
obtain_emiss(out_3st_cont_sim, level = 'subject')
#> $`observation 1`
#> $`observation 1`$`Subject 1`
#>            Mean     SD
#> State 1  49.972 11.778
#> State 2  97.475 11.797
#> State 3 153.374 10.156
#> 
#> $`observation 1`$`Subject 2`
#>            Mean     SD
#> State 1  51.126 11.778
#> State 2 103.527 11.797
#> State 3 150.888 10.156
#> 
#> $`observation 1`$`Subject 3`
#>            Mean     SD
#> State 1  51.562 11.778
#> State 2  95.136 11.797
#> State 3 150.623 10.156
#> 
#> $`observation 1`$`Subject 4`
#>            Mean     SD
#> State 1  50.944 11.778
#> State 2  96.317 11.797
#> State 3 151.328 10.156
#> 
#> $`observation 1`$`Subject 5`
#>            Mean     SD
#> State 1  48.705 11.778
#> State 2  99.862 11.797
#> State 3 151.165 10.156
#> 
#> $`observation 1`$`Subject 6`
#>            Mean     SD
#> State 1  50.108 11.778
#> State 2  98.278 11.797
#> State 3 151.417 10.156
#> 
#> $`observation 1`$`Subject 7`
#>            Mean     SD
#> State 1  49.707 11.778
#> State 2  97.602 11.797
#> State 3 150.817 10.156
#> 
#> $`observation 1`$`Subject 8`
#>            Mean     SD
#> State 1  49.233 11.778
#> State 2  97.115 11.797
#> State 3 146.999 10.156
#> 
#> $`observation 1`$`Subject 9`
#>            Mean     SD
#> State 1  51.440 11.778
#> State 2 101.777 11.797
#> State 3 153.486 10.156
#> 
#> $`observation 1`$`Subject 10`
#>            Mean     SD
#> State 1  51.043 11.778
#> State 2  97.393 11.797
#> State 3 150.500 10.156
#> 
#> 
#> $`observation 2`
#> $`observation 2`$`Subject 1`
#>           Mean    SD
#> State 1  4.613 2.226
#> State 2  9.498 4.954
#> State 3 19.800 3.389
#> 
#> $`observation 2`$`Subject 2`
#>           Mean    SD
#> State 1  6.095 2.226
#> State 2  9.386 4.954
#> State 3 20.436 3.389
#> 
#> $`observation 2`$`Subject 3`
#>           Mean    SD
#> State 1  5.060 2.226
#> State 2  9.831 4.954
#> State 3 20.267 3.389
#> 
#> $`observation 2`$`Subject 4`
#>           Mean    SD
#> State 1  4.773 2.226
#> State 2  8.539 4.954
#> State 3 20.318 3.389
#> 
#> $`observation 2`$`Subject 5`
#>           Mean    SD
#> State 1  5.197 2.226
#> State 2 10.123 4.954
#> State 3 19.810 3.389
#> 
#> $`observation 2`$`Subject 6`
#>           Mean    SD
#> State 1  4.906 2.226
#> State 2  9.505 4.954
#> State 3 20.515 3.389
#> 
#> $`observation 2`$`Subject 7`
#>           Mean    SD
#> State 1  4.513 2.226
#> State 2  9.122 4.954
#> State 3 19.273 3.389
#> 
#> $`observation 2`$`Subject 8`
#>           Mean    SD
#> State 1  4.313 2.226
#> State 2  7.920 4.954
#> State 3 19.664 3.389
#> 
#> $`observation 2`$`Subject 9`
#>           Mean    SD
#> State 1  5.422 2.226
#> State 2  9.566 4.954
#> State 3 19.553 3.389
#> 
#> $`observation 2`$`Subject 10`
#>           Mean    SD
#> State 1  5.233 2.226
#> State 2  9.592 4.954
#> State 3 20.347 3.389

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
