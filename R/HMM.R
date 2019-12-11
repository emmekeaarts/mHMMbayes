#' Onelevel (fixed) hidden  Markov model using Bayesian estimation
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
#' @param one_s_data A matrix containing the observations to be modelled of one
#'   subject or sequence, where the rows represent the observations over time,
#'   and the column(s) contain the dependent variable(s). Note that the
#'   dependent variables have to be numeric, i.e., they cannot be a (set of)
#'   factor variable(s).
#' @param gamma_prior sfsfd
#' @param emiss_prior sf
#' @inheritParams mHMM
#'
#'
#' @return \code{HMM} returns an object of class \code{HMM}, which has
#'   \code{print} and \code{summary} methods to see the results.
#'   The object contains the following components:
#'   \describe{
#'   \item{\code{PD}}{A .. containing ...}
#' }
#'
#' @seealso
#'
#' @examples
#' ### first simulate data
#' # One subject with each n_t observations:
#' n <- 1
#' n_t <- 300
#'
#' # specifying general model properties:
#' m <- 2
#' q_emiss <- 3
#' n_dep <- 1
#' gamma <- matrix(c(0.8, 0.2,
#'                   0.2, 0.8), ncol = m, byrow = TRUE)
#' emiss_distr <- matrix(c(0.50, 0.45, 0.05,
#'                         0.05, 0.05, 0.90), nrow = m, ncol = q_emiss, byrow = TRUE)
#'
#' # simulate data
#' data1 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss,
#'                   gamma = gamma, emiss_distr = emiss_distr)
#'
#' ### analysis of data
#' # specifying starting values
#' st_gamma <- matrix(c(0.7, 0.3,
#'                      0.3, 0.7), ncol = m, byrow = TRUE)
#' st_emiss_distr <- matrix(c(0.40, 0.40, 0.20,
#'                            0.20, 0.20, 0.60), nrow = m, ncol = q_emiss, byrow = TRUE)
#'
#' # run the model
#' out_2st_sim <- HMM(one_s_data = matrix(data1$obs[,2], ncol = 1),
#'                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                     start_val = list(st_gamma, st_emiss_distr),
#'                     mcmc = list(J = 4000, burn_in = 200), show_progress = TRUE)
#'
#' @export
#'
#'

# inlcude check that id is not given as first column
# extent to multivariate
# add description of gamma_prior and emiss_prior
# add see also section to help file
# add references to help file

HMM <- function(one_s_data, gen, start_val, mcmc, return_path = FALSE, show_progress = TRUE,
                 gamma_prior = NULL, emiss_prior = NULL){
  # Initialize data
  n_dep			   <- 1    # gen$n_dep
  m            <- gen$m
  q_emiss 		 <- gen$q_emiss
  dep_labels   <- colnames(one_s_data[,1:(n_dep)])
  J 			     <- mcmc$J
  burn_in			 <- mcmc$burn_in

  # Initalize priors
  if(is.null(gamma_prior)){
    gamma_prior <- matrix(1, ncol = m, nrow = m)
  }
  if(is.null(emiss_prior)){
    emiss_prior <- matrix(1, ncol = q_emiss, nrow = m)
  }

  # Define objects used to store data in mcmc algorithm, not returned
  emiss_counts <- emiss_next <- matrix(0, nrow = m,  ncol = q_emiss)
  gamma_counts <- gamma_next <- matrix(0, nrow = m, ncol = m)

  # Creating a matrix of all for the output of all Bayesian samples (J)
  PD	<- matrix(, nrow = J, ncol = m*q_emiss + m*m + 1)
  colnames(PD) <- c(paste("obs",rep(1:q_emiss, m),"_S",rep(1:m, each = q_emiss), sep = ""), paste("S", rep(1:m, each = m), "to", rep(1:m, m), sep = ""), "LL") #Name giving
  n_t 	<- dim(one_s_data)[1] # Creating a var reflecting the length of the sequence
  sample_path <- matrix(,nrow = n_t, ncol = J) # Creating the matrix wherein all state paths of all Bayesian samples will be stored

  # Put starting values on first line of the outcome matrix
  PD[1, ((sum(m * q_emiss) + 1)) :((sum(m * q_emiss) + m * m))] <- unlist(sapply(start_val, t))[1:(m*m)]
  PD[1, 1:((sum(m * q_emiss)))] <- unlist(sapply(start_val, t))[(m*m + 1): (m*m + sum(m * q_emiss))]

  # start MCMC algorithm
  itime <- proc.time()[3]
  if(show_progress == TRUE){
    cat("Progress of the Bayesian HMM algorithm:", "\n")
    pb <- utils::txtProgressBar(min = 2, max = J, style = 3)
  }
  for (iter in 2 : J){ # For all iterations in all Bayesian samples (J)
    # 1. generate sample path of the Markov chain given s_data, gamma and emission distribution
    emiss 				<- list(matrix(PD[iter-1, 1:(m*q_emiss)], ncol = q_emiss, nrow = m, byrow = TRUE))
    gamma 			  <- matrix(PD[iter-1, (m*q_emiss + 1) : (m*q_emiss + m*m)], byrow = TRUE, ncol = m)
    forward 			<- cat_mult_fw_r_to_cpp(x = one_s_data, m = m, emiss =  emiss, gamma = gamma, n_dep = n_dep, delta = NULL) # obtain forward probabilities
    alpha         <- forward[[1]]
    c             <- max(forward[[2]][, n_t])
    llk           <- c + log(sum(exp(forward[[2]][, n_t] - c)))
    PD[iter, m*q_emiss + m*m + 1] <- llk

    sample_path[n_t, iter] 	<- sample(1:m, 1, prob = alpha[,n_t])
    for (i in (n_t-1):1){
      sample_path[i, iter] <- sample(1:m, 1, prob = (alpha[,i] * gamma[,sample_path[i+1, iter]]))
    }

    # 2. with the sample path, update gamma and conditional distribution
    # update emiss
    for (i in 1:m){
      for (q in 1:q_emiss){
          emiss_counts[i,q] <- sum(sample_path[, iter] == i & one_s_data[,1] == q)
      }
      emiss_next[i,] <- MCMCpack::rdirichlet(1, alpha = (emiss_counts[i,] + emiss_prior[i,]))
    }

    PD[iter,1:(m*q_emiss)] <- as.vector(t(emiss_next))

    # update gamma
    for (i in 1:m){
      for (j in 1:m){
        gamma_counts[i,j] <- sum(sample_path[-n_t, iter] == i & sample_path[-1, iter] == j)
      }
      gamma_next[i,] <- MCMCpack::rdirichlet(1, alpha = (gamma_counts[i,] + gamma_prior[i,]))
    }
    PD[iter,(m*q_emiss+1):(m*q_emiss+m*m)] <- as.vector(t(gamma_next))

    if(show_progress == TRUE){
      utils::setTxtProgressBar(pb, iter)
    }
  }

  # End of function, return output values --------
  if(show_progress == TRUE){
    close(pb)
  }
  ctime = proc.time()[3]
  message(paste("Total time elapsed (hh:mm:ss):", hms(ctime-itime)))

  if(return_path == TRUE){
    out <- list(input = list(m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,
                             burn_in = burn_in, n_t = n_t, dep_labels = dep_labels),
                PD = PD, sample_path= sample_path)
  } else {
    out <- list(input = list(m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,
                             burn_in = burn_in, n_t = n_t, dep_labels = dep_labels),
                PD = PD)
  }
  class(out) <- append(class(out), "HMM")
  return(out)
}


