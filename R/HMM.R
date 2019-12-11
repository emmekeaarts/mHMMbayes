#' Onelevel (fixed) hidden  Markov model using Bayesian estimation
#'
#' \code{HMM} fits a onelevel (i.e., fixed parameter) hidden Markov model (HMM)
#' to intense longitudinal data with categorical observations using Bayesian
#' estimation, and creates an object of class HMM.
#'
#' For a short description of
#' the package see \link{mHMMbayes}. See \code{vignette("tutorial-mhmm")} for an
#' introduction to multilevel hidden Markov models and the package, and see
#' \code{vignette("estimation-mhmm")} for an overview of the used estimation
#' algorithms.
#'
#'
#' @param s_data A matrix containing the observations to be modelled of one
#'   subject or sequence, where the rows represent the observations over time.
#'   In \code{s_data}, the first column indicates subject id number, which
#'   should be identitical over the rows for \code{HMM()}. The subsequent
#'   columns contain the dependent variable(s). Hence, the number of columns are
#'   equal to the number of dependent variables (\code{n_dep}) + 1. Note that
#'   the dependent variables have to be numeric, i.e., they cannot be a (set of)
#'   factor variable(s).
#' @param gamma_prior A \code{m} by \code{m} matrix containing the parameters of
#'   the Dirichlet prior distribution on the transition probability matrix. When
#'   left unspecified, defaults to an uninformative prior (that is, all entries
#'   of the matrix equal 1).
#' @param emiss_prior A \code{m} by \code{q_emiss} matrix containing the parameters of
#'   the Dirichlet prior distribution on the categrical emission distribution. When
#'   left unspecified, defaults to an uninformative prior (that is, all entries
#'   of the matrix equal 1).
#' @inheritParams mHMM
#'
#'
#' @return \code{HMM} returns an object of class \code{HMM}, which has
#'   \code{print} and \code{summary} methods to see the results.
#'   The object contains the following components:
#'   \describe{
#'   \item{\code{PD}}{A matrix containing the posterior distribution of the
#'   parameter estimates and the log likelihood over the iterations of the MCMC
#'   sampler. The iterations of the sampler are contained in the rows, and the
#'   columns contain the estimates of subsequently the emission probabilities,
#'   the transition probabilities and the log likelihood.}
#'   \item{\code{input}}{Overview of used input specifications: the number of
#'   states \code{m}, the number of used dependent variables \code{n_dep}, the
#'   number of output categories for each of the dependent variables
#'   \code{q_emiss}, the number of iterations \code{J} and the specified burn in
#'   period \code{burn_in} of the MCMC sampler, the observation length
#'   \code{n_t}, and the column names of the dependent variables
#'   \code{dep_labels}.}
#'   \item{\code{sample_path}}{A matrix containing the sampled hidden state
#'   sequence over the MCMC sampler. The time points of the dataset are
#'   contained in the rows, and the sampled paths over the iterations are
#'   contained in the columns. Only returned if \code{return_path = TRUE}. }
#' }
#'
#' @seealso \code{\link{sim_mHMM}} for simulating (multilevel) hidden Markov data,
#'   \code{\link{vit_mHMM}} for obtaining the most likely hidden state sequence
#'   using the Viterbi algorithm, \code{\link{obtain_gamma}}
#'   and \code{\link{obtain_emiss}} for obtaining the transition or emission
#'   distribution probabilities of a fitted model.
#'
#' @references
#' \insertRef{rabiner1989}{mHMMbayes}
#'
#' \insertRef{scott2002}{mHMMbayes}
#'
#' \insertRef{zucchini2017}{mHMMbayes}
#'
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
#' out_2st_sim <- HMM(one_s_data = matrix(data1$obs, ncol = 1),
#'                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                     start_val = list(st_gamma, st_emiss_distr),
#'                     mcmc = list(J = 4000, burn_in = 200), show_progress = TRUE)
#'
#' @export
#'
#'
# help file uitbreiden, korte beschrijving Bayes methode (bij bayes combineren we prior en likelihood,
# gebruiken hier Dirichlet prior, waardes van Dirichlet distributie kunnen geinterpreteerd worden
# als .., forward probabilities, backward sampling).
# create plotting function plot.HMM() for posterior densities for HMM as well (see plot.mHMM())
# include more checks (i.e., dimensions of priors if given, etc)
# gamma_next and emiss_next are obsolete, delete.
# vit_mHMM compatible maken
# obtain_gamma en obtain_emiss compatible maken

HMM <- function(s_data, gen, start_val, mcmc, return_path = FALSE, show_progress = TRUE,
                 gamma_prior = NULL, emiss_prior = NULL){
  # Initialize data
  if(length(unique(s_data[,1])) != 1 ){
    stop("The first column in the dataset provided should contain the subject / sequence ID, which should be identical over the rows for use of HMM(). For analyzing multiple sequences at once, please use mHMM() instead.")
  }
  n_dep			   <- gen$n_dep
  if(dim(s_data)[2] != n_dep + 1 ){
    stop("The number of columns in the dataset provided should equal the number of dependent variables n_dep + 1.")
  }
  one_s_data   <- s_data[,-1]
  m            <- gen$m
  q_emiss 		 <- gen$q_emiss
  start        <- c(0, q_emiss * m)
  dep_labels   <- colnames(one_s_data[,1:(n_dep)])
  n_t 	       <- dim(one_s_data)[1] # Creating a var reflecting the length of the sequence
  J 			     <- mcmc$J
  burn_in			 <- mcmc$burn_in

  # Initalize priors
  if(is.null(gamma_prior)){
    gamma_prior <- matrix(1, ncol = m, nrow = m)
  }
  if(is.null(emiss_prior)){
    emiss_prior <- lapply(q_emiss, dif_matrix, rows = m, data = 1)
  }

  # Define objects used to store data in mcmc algorithm, not returned
  emiss        <- rep(list(NULL), n_dep)
  emiss_counts <- emiss_next <- lapply(q_emiss, dif_matrix, rows = m, data = 0)
  gamma_counts <- gamma_next <- matrix(0, nrow = m, ncol = m)

  # Creating a matrix of all for the output of all Bayesian samples (J)
  PD	<- matrix(, nrow = J, ncol = sum(m*q_emiss) + m*m + 1)
  PD_emiss_names   <- paste("q", 1, "_emiss", rep(1:q_emiss[1], m), "_S", rep(1:m, each = q_emiss[1]), sep = "")
  if(n_dep > 1){
    for(q in 2:n_dep){
      PD_emiss_names <- c(PD_emiss_names, paste("q", q, "_emiss", rep(1:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q]), sep = ""))
    }
  }
  colnames(PD) 	<- c(PD_emiss_names, paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = ""), "LL")
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
    for(q in 1:n_dep){
      emiss[[q]] <- matrix(PD[iter-1,(sum(start[1:q]) + 1):(sum(start[1:q]) + (m * q_emiss[q]))], byrow = TRUE, ncol = q_emiss[q], nrow = m)
    }
    gamma 			  <- matrix(PD[iter-1, (sum(m*q_emiss) + 1) : (sum(m*q_emiss) + m*m)], byrow = TRUE, ncol = m)
    forward 			<- cat_mult_fw_r_to_cpp(x = one_s_data, m = m, emiss =  emiss, gamma = gamma, n_dep = n_dep, delta = NULL) # obtain forward probabilities
    alpha         <- forward[[1]]
    c             <- max(forward[[2]][, n_t])
    llk           <- c + log(sum(exp(forward[[2]][, n_t] - c)))
    PD[iter, sum(m * q_emiss) + m * m + 1] <- llk

    sample_path[n_t, iter] 	<- sample(1:m, 1, prob = alpha[,n_t])
    for (i in (n_t-1):1){
      sample_path[i, iter] <- sample(1:m, 1, prob = (alpha[,i] * gamma[,sample_path[i+1, iter]]))
    }

    # 2. with the sample path, update gamma and conditional distribution
    # update emiss
    for (q in 1:n_dep){
      for (i in 1:m){
        for (k in 1:q_emiss[q]){
          emiss_counts[[q]][i,k] <- sum(sample_path[, iter] == i & one_s_data[,q] == k)
        }
        emiss_next[[q]][i,] <- PD[iter, (sum(start[1:q]) + 1 + (i - 1) * q_emiss[q]):(sum(start[1:q]) + (i - 1) * q_emiss[q] + q_emiss[q])] <-
          MCMCpack::rdirichlet(1, alpha = (emiss_counts[[q]][i,] + emiss_prior[[q]][i,]))
      }
    }

    # update gamma
    for (i in 1:m){
      for (j in 1:m){
        gamma_counts[i,j] <- sum(sample_path[-n_t, iter] == i & sample_path[-1, iter] == j)
      }
      gamma_next[i,] <- MCMCpack::rdirichlet(1, alpha = (gamma_counts[i,] + gamma_prior[i,]))
    }
    PD[iter,(sum(m * q_emiss) + 1):(sum(m*q_emiss) + m * m)] <- as.vector(t(gamma_next))

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


