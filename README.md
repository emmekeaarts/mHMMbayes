
<!-- README.md is generated from README.Rmd. Please edit that file -->
mHMMbayes
=========

With the package mHMMbayes you can fit multilevel hidden Markov models. The multilevel hidden Markov model (HMM) is a generalization of the well-known hidden Markov model, tailored to accomodate (intense) longitudinal data of multiple individuals simultaneously. Using a multilevel framework, we allow for heterogeneity in the model parameters (transition probability matrix and conditional distribution), while estimating one overall HMM. The model has a great potential of application in many fields, such as the social sciences and medicine. The model can be fitted on multivariate data with a catagorical distribution, and include individual level covariates (allowing for e.g., group comparisons on model parameters). Parameters are estimated using Bayesian estimation utilizing the forward-backward recursion within a hybrid Metropolis within Gibbs sampler. The package also includes a function to simulate data and a function to obtain the most likely hidden state sequence for each individual using the Viterbi algorithm.

NOTE: this is a beta version of the package. The package is still (heavily) under development, new functionalities will be added and functionalities will be altered in the near future.

Please do not hesitate to contact me if you have any questions regarding the package.

Installation
------------

You can install mHMMbayes from github with:

``` r
# install.packages("devtools")
devtools::install_github("emmekeaarts/mHMMbayes")
```

Usage
-----

This is a basic example which shows you how to run the model using example data included with the package, and how to simulate data:

``` r
library(mHMMbayes)
# specifying general model properties
m <- 2
n_dep <- 4
q_emis <- c(3, 2, 3, 2)

# specifying starting values
start.EM <- list(matrix(c(0.9, 0.05, 0.05, 0.05, 0.05, 0.9), byrow = TRUE,
                         nrow = m, ncol = q_emis[1]), # vocalizing patient
                  matrix(c(0.9, 0.1, 0.9, 0.1), byrow = TRUE, nrow = m,
                         ncol = q_emis[2]), # looking patient
                  matrix(c(0.05, 0.05, 0.9, 0.9, 0.05, 0.05), byrow = TRUE,
                         nrow = m, ncol = q_emis[3]), # vocalizing therapist
                  matrix(c(0.9, 0.1, 0.9, 0.1), byrow = TRUE, nrow = m,
                         ncol = q_emis[4])) # looking therapist
 start.TM <- diag(.8, m)
 start.TM[lower.tri(start.TM) | upper.tri(start.TM)] <- .2

 # run a model without covariates
 set.seed(23245)
 out <- mHMM_mnl(s_data = nonverbal, gen = list(m = m, n_dep = n_dep,
                 q_emis = q_emis), start_val = c(as.vector(t(start.EM[[1]])),
                 as.vector(t(start.EM[[2]])), as.vector(t(start.EM[[3]])),
                 as.vector(t(start.EM[[4]])), as.vector(t(start.TM))),
                 mcmc = list(J = 11, burn_in = 5))
#> [1] 10
#> [1] "total time elapsed in minutes 0.37"

 # including covariates. Only the emission distribution for each of the 4
 # dependent variables is predicted using standardized CDI change.
 n_subj <- 10
 xx <- rep(list(matrix(1, ncol = 1, nrow = n_subj)), (n_dep + 1))
 for(i in 2:(n_dep + 1)){
   xx[[i]] <- cbind(xx[[i]], nonverbal_cov$std_CDI_change)
 }
 set.seed(34109)
 out2 <- mHMM_mnl(s_data = nonverbal, xx = xx, gen = list(m = m, n_dep = n_dep,
                 q_emis = q_emis), start_val = c(as.vector(t(start.EM[[1]])),
                 as.vector(t(start.EM[[2]])), as.vector(t(start.EM[[3]])),
                 as.vector(t(start.EM[[4]])), as.vector(t(start.TM))),
                 mcmc = list(J = 11, burn_in = 5))
#> [1] 10
#> [1] "total time elapsed in minutes 0.36"

 
 ### Simulating data
 # simulating data for 10 subjects with each 100 observations
 T <- 100
 n <- 10
 m <- 3
 pr <- 4
 gamma <- matrix(c(0.8, 0.1, 0.1,
                   0.2, 0.7, 0.1,
                   0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
 emiss_distr <- matrix(c(0.5, 0.5, 0.0, 0.0,
                         0.1, 0.1, 0.8, 0.0,
                         0.0, 0.0, 0.1, 0.9), nrow = m, ncol = pr, byrow = TRUE)
 set.seed(1253)
 data1 <- sim_mHMM(T = T, n = n, m = m, pr = pr, gamma = gamma, emiss_distr = emiss_distr,
                   var_gamma = 1, var_emiss = 1)


 # simulating subject specific transition probability matrices and emission distributions only
 T <- 0
 n <- 5
 m <- 3
 pr <- 4
 gamma <- matrix(c(0.8, 0.1, 0.1,
                   0.2, 0.7, 0.1,
                   0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
 emiss_distr <- matrix(c(0.5, 0.5, 0.0, 0.0,
                         0.1, 0.1, 0.8, 0.0,
                         0.0, 0.0, 0.1, 0.9), nrow = m, ncol = pr, byrow = TRUE)
 set.seed(549801)
 data2 <- sim_mHMM(T = T, n = n, m = m, pr = pr, gamma = gamma, emiss_distr = emiss_distr,
                   var_gamma = 1, var_emiss = 1)
 data2
#> $subject_gamma
#> $subject_gamma[[1]]
#>        [,1]   [,2]   [,3]
#> [1,] 0.6302 0.2849 0.0849
#> [2,] 0.1817 0.7714 0.0469
#> [3,] 0.2164 0.1738 0.6098
#> 
#> $subject_gamma[[2]]
#>        [,1]   [,2]   [,3]
#> [1,] 0.7819 0.1235 0.0945
#> [2,] 0.0747 0.9015 0.0238
#> [3,] 0.3285 0.4705 0.2011
#> 
#> $subject_gamma[[3]]
#>        [,1]   [,2]   [,3]
#> [1,] 0.6228 0.1443 0.2329
#> [2,] 0.5242 0.3106 0.1652
#> [3,] 0.5215 0.1167 0.3618
#> 
#> $subject_gamma[[4]]
#>        [,1]   [,2]   [,3]
#> [1,] 0.5726 0.1054 0.3220
#> [2,] 0.1751 0.5438 0.2811
#> [3,] 0.2109 0.1686 0.6204
#> 
#> $subject_gamma[[5]]
#>        [,1]   [,2]   [,3]
#> [1,] 0.8227 0.1212 0.0561
#> [2,] 0.2029 0.5990 0.1982
#> [3,] 0.0902 0.3200 0.5898
#> 
#> 
#> $subject_emmis
#> $subject_emmis[[1]]
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.4916 0.5083 0.0001 0.0000
#> [2,] 0.0572 0.1810 0.7618 0.0000
#> [3,] 0.0000 0.0000 0.0629 0.9371
#> 
#> $subject_emmis[[2]]
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.1451 0.8549 0.0000 0.0000
#> [2,] 0.0510 0.1518 0.7972 0.0000
#> [3,] 0.0000 0.0000 0.0514 0.9486
#> 
#> $subject_emmis[[3]]
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.2158 0.7842 0.0000 0.0000
#> [2,] 0.1865 0.1831 0.6304 0.0000
#> [3,] 0.0001 0.0002 0.5793 0.4204
#> 
#> $subject_emmis[[4]]
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.3268 0.6732 0.0000 0.0000
#> [2,] 0.0501 0.0867 0.8631 0.0000
#> [3,] 0.0000 0.0000 0.1194 0.8806
#> 
#> $subject_emmis[[5]]
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.7268 0.2731 0.0000 0.0001
#> [2,] 0.1303 0.2549 0.6148 0.0000
#> [3,] 0.0000 0.0000 0.1085 0.8915

 set.seed(10893)
 data3 <- sim_mHMM(T = T, n = n, m = m, pr = pr, gamma = gamma, emiss_distr = emiss_distr,
                   var_gamma = .5, var_emiss = .5)
 data3
#> $subject_gamma
#> $subject_gamma[[1]]
#>        [,1]   [,2]   [,3]
#> [1,] 0.7958 0.0461 0.1581
#> [2,] 0.3042 0.4663 0.2295
#> [3,] 0.1396 0.6520 0.2084
#> 
#> $subject_gamma[[2]]
#>        [,1]   [,2]   [,3]
#> [1,] 0.6221 0.0834 0.2946
#> [2,] 0.1143 0.8430 0.0427
#> [3,] 0.1414 0.3805 0.4780
#> 
#> $subject_gamma[[3]]
#>        [,1]   [,2]  [,3]
#> [1,] 0.7416 0.0554 0.203
#> [2,] 0.1403 0.7937 0.066
#> [3,] 0.1915 0.0975 0.711
#> 
#> $subject_gamma[[4]]
#>        [,1]   [,2]   [,3]
#> [1,] 0.6333 0.0932 0.2736
#> [2,] 0.1127 0.6909 0.1964
#> [3,] 0.1058 0.3872 0.5070
#> 
#> $subject_gamma[[5]]
#>        [,1]   [,2]   [,3]
#> [1,] 0.7610 0.1833 0.0557
#> [2,] 0.0781 0.8920 0.0300
#> [3,] 0.2269 0.1116 0.6615
#> 
#> 
#> $subject_emmis
#> $subject_emmis[[1]]
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.6359 0.3641 0.0000 0.0000
#> [2,] 0.2302 0.2625 0.5073 0.0000
#> [3,] 0.0000 0.0000 0.0326 0.9674
#> 
#> $subject_emmis[[2]]
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.6471 0.3528 0.0000 0.0000
#> [2,] 0.1977 0.2782 0.5240 0.0001
#> [3,] 0.0000 0.0000 0.2446 0.7554
#> 
#> $subject_emmis[[3]]
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.7053 0.2946 0.0000 0.0000
#> [2,] 0.1433 0.0626 0.7940 0.0001
#> [3,] 0.0000 0.0000 0.0274 0.9726
#> 
#> $subject_emmis[[4]]
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.8119 0.1880 0.0000 0.0000
#> [2,] 0.0704 0.0669 0.8627 0.0000
#> [3,] 0.0000 0.0000 0.1557 0.8443
#> 
#> $subject_emmis[[5]]
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.5968 0.4032 0.0000 0.0000
#> [2,] 0.1023 0.1158 0.7819 0.0000
#> [3,] 0.0000 0.0000 0.1358 0.8642
```
