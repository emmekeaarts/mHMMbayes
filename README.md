
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
The model can be fitted on multivariate data with a categorical
distribution, and include individual level covariates (allowing for
e.g., group comparisons on model parameters). Parameters are estimated
using Bayesian estimation utilizing the forward-backward recursion
within a hybrid Metropolis within Gibbs sampler. The package also
includes various options for model visualization, a function to simulate
data and a function to obtain the most likely hidden state sequence for
each individual using the Viterbi algorithm.

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
out_3st_vary_dep <- mHMM_vary(s_data = data_vary$obs,
                    gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                    data_distr = data_distr,
                    start_val = c(list(gamma), emiss_distr),
                    emiss_cont_hyp_prior = hyp_pr,
                    mcmc = list(J = 11, burn_in = 5))
#> Progress of the Bayesian mHMM algorithm: 
#>   |                                                                              |                                                                      |   0%  |                                                                              |========                                                              |  11%  |                                                                              |================                                                      |  22%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================                                       |  44%  |                                                                              |=======================================                               |  56%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================                |  78%  |                                                                              |==============================================================        |  89%  |                                                                              |======================================================================| 100%
#> Total time elapsed (hh:mm:ss): 00:00:00

out_3st_vary_dep
#> Number of subjects: 10 
#> 
#> 11 iterations used in the MCMC algorithm with a burn in of 5 
#> Average Log likelihood over all subjects: -402.8316 
#> Average AIC over all subjects: 853.6631 
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
#> From state 1      0.516      0.382      0.097
#> From state 2      0.212      0.647      0.121
#> From state 3      0.165      0.369      0.462
#> 
#>  
#> Emission distribution for each of the dependent variables 
#>  (at the group level): 
#>  
#> $`observation 1`
#>         Category 1 Category 2 Category 3 Category 4
#> State 1      0.744      0.177      0.060      0.026
#> State 2      0.104      0.765      0.080      0.026
#> State 3      0.062      0.270      0.377      0.273
#> 
#> $`observation 2`
#>          Mean Variance
#> State 1 5.866    0.560
#> State 2 5.283    0.632
#> State 3 2.895    0.905
#> 
#> $`observation 3`
#>          Mean Variance
#> State 1 3.955    1.264
#> State 2 5.310    1.229
#> State 3 6.152    1.398
```
