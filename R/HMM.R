#' Onelevel hidden  Markov model using Bayesian estimation
#'
#' \code{HMM} fits a onelevel (i.e., fixed parameter) hidden Markov model (HMM)
#' to intense longitudinal data with categorical observations using Bayesian
#' estimation, and creates an object of class HMM. For a short description of
#' the package see \link{mHMMbayes}. See \code{vignette("tutorial-mhmm")} for an
#' introduction to multilevel hidden Markov models and the package, and see
#' \code{vignette("estimation-mhmm")} for an overview of the used estimation
#' algorithms.
#'
#'
#' @param s_data A matrix containing the observations to be modelled, where the
#'   rows represent the observations over time. In \code{s_data}, the first
#'   column indicates subject id number. Hence, the id number is repeated over
#'   rows equal to the number of observations for that subject. The subsequent
#'   columns contain the dependent variable(s). Note that the dependent
#'   variables have to be numeric, i.e., they cannot be a (set of) factor
#'   variable(s). The total number of rows are equal to the sum over the number
#'   of observations of each subject, and the number of columns are equal to the
#'   number of dependent variables (\code{n_dep}) + 1. The number of
#'   observations can vary over subjects.
#' @param xx sfsfd
#'
#'
#' @return \code{HMM} returns an object of class \code{HMM}, which has
#'   \code{print} and \code{summary} methods to see the results.
#'   The object contains the following components:
#'   \describe{
#'   \item{\code{xx}}{A .. containing ...}
#' }
#'
#' @seealso
#'
#' @examples
#'
#' @export
#'
#'

# Incorporate show_progress
# inlcude check that id is not given as first column
# extent to multivariate
# replace alpha.cond.dist0.sc, alpha.cond.dist0 and alpha.gamma0 with emiss_prior and gamma_prior

HMM <- function(s_data, gen, start_val, mcmc, return_path = FALSE, show_progress = TRUE,
                 gamma_prior = NULL, emiss_prior = NULL){
  n_dep			   <- 1    # gen$n_dep
  m            <- gen$m
  q_emiss 		 <- gen$q_emiss
  PD	<- matrix(, nrow = J, ncol = m*q_emiss + m*m) # Creating a matrix of all Bayesian samples (J) by transition and emission probs.
  colnames(PD) <- c(paste("obs",rep(1:q_emiss, m),"_S",rep(1:m, each = q_emiss), sep = ""), paste("S", rep(1:m, each = m), "to", rep(1:m, m), sep = "")) #Name giving
  n 	<- dim(s_data)[1] # Creating a var reflecting the length of the sequence
  sample_path <- matrix(,nrow = n, ncol = J) # Creating the matrix wherein all state paths of all Bayesian samples will be stored

  alpha.cond.dist0.sc <- alpha.cond.dist0 #Copying the initial state distribution with only 1's.

    # Put starting values on first line of the outcome matrix
  PD[1, ((sum(m * q_emiss) + 1)) :((sum(m * q_emiss) + m * m))] <- unlist(sapply(start_val, t))[1:(m*m)]
  PD[1, 1:((sum(m * q_emiss)))] <- unlist(sapply(start_val, t))[(m*m + 1): (m*m + sum(m * q_emiss))]

  for (iter in 2 : J){ # For all iterations in all Bayesian samples (J)
    # 1. generate sample path of the Markov chain given s_data, gamma and emission distribution
    emiss 				<- matrix(PD[iter-1, 1:(m*q_emiss)], ncol = q_emiss, nrow = m, byrow = TRUE)
    gamma 			  <- matrix(PD[iter-1, (m*q_emiss + 1) : (m*q_emiss + m*m)], byrow = TRUE, ncol = m)
    forward 			<- cat_mult_fw_r_to_cpp(x = s_data, m = m, emiss =  emiss, gamma = gamma, delta = NULL) # obtain forward probabilities
    alpha         <- forward[[1]]
    c             <- max(forward[[2]][, subj_data[[s]]$n])
    # INCORPORATE IN OUTPUT!
    llk           <- c + log(sum(exp(forward[[2]][, subj_data[[s]]$n] - c)))

    sample_path[n, iter] 	<- sample(1:m, 1, prob = alpha[,n])
    for (i in (n-1):1){
      sample_path[i, iter] <- sample(1:m, 1, prob = (alpha[,i] * gamma[,sample_path[i+1, iter]]))
    }

    # 2. with the sample path, update gamma and conditional distribution
    # update cond.dist
    cond.dist.counts <- matrix(0,m,q_emiss)
    for (i in 2:n){
      for (j in 1:m){
        for (k in 1:q_emiss){
          if (sample_path[i, iter] == j & s_data[i] == k) {
            cond.dist.counts[j,k] <- cond.dist.counts[j,k] + 1
          }
        }
      }
    }
    cond.dist.next <- matrix(, m, q_emiss)
    for (i in 1:m){
      cond.dist.next[i,] <- rdirichlet(1, alpha = (cond.dist.counts[i,] + alpha.cond.dist0.sc[i,]))
    }
    PD[iter,1:(m*q_emiss)] <- as.vector(t(cond.dist.next))

    # update gamma COUNTS DING MAKEN
    gamma.counts <- matrix(0, m, m)
    for (i in 2:n){
      for (j in 1:m){
        for (k in 1:m){
          if (sample_path[i-1, iter] == j & sample_path[i, iter] == k) {
            gamma.counts[j,k] <- gamma.counts[j,k] + 1
          }
        }
      }
    }
    gamma.next <- matrix(, m, m)
    for (i in 1:m){
      gamma.next[i,] <- rdirichlet(1, alpha = (gamma.counts[i,] + alpha.gamma0))
    }
    PD[iter,(m*q_emiss+1):(m*q_emiss+m*m)] <- as.vector(t(gamma.next))
  }
  return(list(PD = PD, sample_path= sample_path))

}


