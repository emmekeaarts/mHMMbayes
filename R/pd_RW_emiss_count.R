#' Proposal distribution settings RW Metropolis sampler for mHMM
#' Poisson-lognormal emission distribution(s)
#'
#' \code{pd_RW_emiss_count} provides a framework to manually specify the
#' settings of the proposal distribution of the random walk (RW) Metropolis
#' sampler of the emission distribution(s) of the multilevel hidden Markov
#' model, and creates on object of the class \code{mHMM_pdRW_emiss}. The RW
#' metropolis sampler is used for sampling the subject level parameter estimates
#' relating to the emission distributions of the dependent variables \code{k},
#' that is, the Poisson parameters lambda.
#'
#'
#' When no manual values for the settings of the proposal distribution of the
#' random walk (RW) Metropolis sampler are specified at all (that is, the
#' function \code{pd_RW_emiss_count} is not used), \code{emiss_scalar} set
#' to 2.38, and \code{emiss_w} set to 0.1. See the section \emph{Scaling the
#' proposal distribution of the RW Metropolis sampler} in
#' \code{vignette("estimation-mhmm")} for details.
#'
#' Within the function \code{mHMM}, the acceptance rate of the RW metropolis
#' sampler relating to the emission distribution(s) can be tracked using the
#' output parameter \code{emiss_naccept}. An acceptance rate of about 45\% is
#' considered optimal when a single parameter is being updated (Gelman,
#' Carlin, Stern & Rubin, 2014).
#'
#' @inheritParams mHMM
#' @param emiss_scalar A list containing \code{n_dep} elements
#'  corresponding to each of the dependent variables, where each element is a
#'  numeric vector with length 1 denoting the scale factor \code{s}. That is,
#'  the scale of the proposal distribution is composed of a covariance matrix
#'  Sigma, which is then tuned by multiplying it by a scaling factor \code{s}^2.
#' @param emiss_w A list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables, where each element is a numeric vector
#'  with length 1 denoting the weight for the overall log likelihood (i.e., log
#'  likelihood based on the pooled data over all subjects) in the fractional
#'  likelihood.
#'
#' @return \code{pd_RW_emiss_count} returns an object of class
#'   \code{mHMM_pdRW_emiss}, containing settings of the proposal distribution of
#'   the random walk (RW) Metropolis sampler on the categorical emission
#'   distribution(s) of the multilevel hidden Markov model. The object is
#'   specifically created and formatted for use by the function \code{mHMM}, and
#'   checked for correct input dimensions. The object contains the following
#'   components:
#'    \describe{
#'   \item{\code{gen}}{A list containing the elements \code{m} and
#'   \code{n_dep}, used for checking equivalent general model properties
#'   specified under \code{pd_RW_emiss_count} and \code{mHMM}.}
#'   \item{\code{emiss_scalar}}{A list containing \code{n_dep} elements denoting
#'   the scale factor \code{s} of the proposal distribution.}
#'   \item{\code{emiss_w}}{A list containing \code{n_dep} elements denoting
#'   denoting the weight for the overall log likelihood in the fractional
#'   likelihood.}
#'   }
#'
#' @references
#' \insertRef{gelman2014}{mHMMbayes}
#'
#' \insertRef{rossi2012}{mHMMbayes}
#'
#'
#' @examples
#' ###### Example using package simulated data
#' # specifying general model properties:
#' m <- 3
#' n_dep <- 3
#'
#' # specifying manual values for RW metropolis sampler on emission distributions
#' emiss_mle0 <- list(matrix(c( 2,  0,
#'                                 -2, -2,
#'                                  0, -1), byrow = TRUE, nrow = m, ncol = q_emiss[1] - 1),
#'                        matrix(c( 2,
#'                                  2,
#'                                  2), byrow = TRUE, nrow = m, ncol = q_emiss[2] - 1),
#'                        matrix(c(-2, -2,
#'                                  2,  0,
#'                                  0, -1), byrow = TRUE, nrow = m, ncol = q_emiss[3] - 1),
#'                        matrix(c( 2,
#'                                  2,
#'                                  2), byrow = TRUE, nrow = m, ncol = q_emiss[4] - 1))
#' emiss_scalar <- list(c(2), c(3))
#' emiss_w <- rep(list(c(0.2)), n_dep)
#' manual_emiss_sampler <- pd_RW_emiss_count(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                                         emiss_mle0 = emiss_mle0,
#'                                         emiss_scalar = emiss_scalar,
#'                                         emiss_w = emiss_w)
#'
#' # specifying starting values
#' start_TM <- diag(.7, m)
#' start_TM[lower.tri(start_TM) | upper.tri(start_TM)] <- .1
#' start_EM <- list(matrix(c(0.05, 0.90, 0.05,
#'                           0.90, 0.05, 0.05,
#'                           0.55, 0.45, 0.05), byrow = TRUE,
#'                         nrow = m, ncol = q_emiss[1]), # vocalizing patient
#'                  matrix(c(0.1, 0.9,
#'                           0.1, 0.9,
#'                           0.1, 0.9), byrow = TRUE, nrow = m,
#'                         ncol = q_emiss[2]), # looking patient
#'                  matrix(c(0.90, 0.05, 0.05,
#'                           0.05, 0.90, 0.05,
#'                           0.55, 0.45, 0.05), byrow = TRUE,
#'                         nrow = m, ncol = q_emiss[3]), # vocalizing therapist
#'                  matrix(c(0.1, 0.9,
#'                           0.1, 0.9,
#'                           0.1, 0.9), byrow = TRUE, nrow = m,
#'                         ncol = q_emiss[4])) # looking therapist
#'
#' # Note that for reasons of running time, J is set at a ridiculous low value.
#' # One would typically use a number of iterations J of at least 1000,
#' # and a burn_in of 200.
#' \donttest{
#' out_3st_RWemiss <- mHMM(s_data = nonverbal,
#'                          gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                          start_val = c(list(start_TM), start_EM),
#'                          emiss_sampler = manual_emiss_sampler,
#'                          mcmc = list(J = 11, burn_in = 5))
#'
#' out_3st_RWemiss
#' summary(out_3st_RWemiss)
#'
#' # checking acceptance rate (for illustrative purposes, in the example,
#' # J is too low for getting a fair indication)
#' div_J <- function(x, J) x / J
#' J_it <- 11 - 1 # accept/reject starts at iteration 2 of MCMC algorithm
#' RW_emiss_accept <- sapply(out_3st_RWemiss$emiss_naccept, div_J, J_it, simplify = FALSE)
#'
#' # average acceptance rate over all subjects per parameter
#' # rows represent each of the n_dep dependent variables, columns represent the m states
#' t(sapply(RW_emiss_accept, apply, MARGIN = 2, mean))
#' }
#' \dontshow{
#' out_3st_RWemiss <- mHMM(s_data = nonverbal,
#'                          gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                          start_val = c(list(start_TM), start_EM),
#'                          emiss_sampler = manual_emiss_sampler,
#'                          mcmc = list(J = 5, burn_in = 3))
#' }
#' n_t     <- 200     # Number of observations on the dependent variable
#' m       <- 3        # Number of hidden states
#' n_dep   <- 3        # Number of dependent variables
#' n_subj  <- 30        # Number of subjects
#'
#' gamma   <- matrix(c(0.9, 0.05, 0.05,
#'                     0.2, 0.7, 0.1,
#'                     0.2,0.3, 0.5), ncol = m, byrow = TRUE)
#'
#' emiss_distr <- list(matrix(c(log(20),
#'                              log(10),
#'                              log(5)), nrow = m, byrow = TRUE),
#'                     matrix(c(log(15),
#'                              log(2),
#'                              log(5)), nrow = m, byrow = TRUE),
#'                     matrix(c(log(50),
#'                              log(3),
#'                              log(20)), nrow = m, byrow = TRUE))
#'
#' # Simulate count data
#' data_count <- sim_mHMM(n_t = n_t, n = n_subj,
#'                        data_distr = "count", gen = list(m = m, n_dep = n_dep),
#'                        gamma = gamma, emiss_distr = emiss_distr,
#'                        var_gamma = 0.1, var_emiss = rep(0.01, n_dep), return_ind_par = TRUE)
#'
#' # Transition probabilities
#' start_gamma <- diag(0.8, m)
#' start_gamma[lower.tri(start_gamma) | upper.tri(start_gamma)] <- (1 - diag(start_gamma)) / (m - 1)
#'
#' # Emission distribution
#' start_emiss <- list(matrix(c(20,10, 5), nrow = m, byrow = TRUE),
#'                     matrix(c(15, 2, 5), nrow = m, byrow = TRUE),
#'                     matrix(c(50, 3,20), nrow = m, byrow = TRUE))
#'
#' # Specify hyper-prior for the count emission distribution
#' manual_prior_emiss <- prior_emiss_count(
#'   gen = list(m = m, n_dep = n_dep),
#'   emiss_mu0 = list(matrix(log(c(20, 10, 5)), byrow = TRUE, ncol = m),
#'                    matrix(log(c(15, 2, 5)), byrow = TRUE, ncol = m),
#'                    matrix(log(c(50, 3, 20)), byrow = TRUE, ncol = m)),
#'   emiss_K0  = rep(list(0.1),n_dep),
#'   emiss_nu  = rep(list(0.1),n_dep),
#'   emiss_V   = rep(list(rep(10, m)),n_dep)
#' )
#'
#' # Specify the desired values for the sampler
#' manual_emiss_sampler <- pd_RW_emiss_count(gen = list(m = m, n_dep = n_dep),
#'                                         emiss_scalar = rep(list(2.38),n_dep),
#'                                         emiss_w = rep(list(0.1), n_dep))
#'
#' # Run model
#' # Note that for reasons of running time, J is set at a ridiculous low value.
#' # One would typically use a number of iterations J of at least 1000,
#' # and a burn_in of 200.
#' out_3st_count_sim <- mHMM(s_data = data_count$obs,
#'                           data_distr = 'count',
#'                           gen = list(m = m, n_dep = n_dep),
#'                           start_val = c(list(start_gamma), start_emiss),
#'                           emiss_hyp_prior = manual_prior_emiss,
#'                           emiss_sampler = manual_emiss_sampler,
#'                           mcmc = list(J = 11, burn_in = 5),
#'                           show_progress = TRUE)
#'
#' # Examine acceptance rates over dependent variable, individual, and states:
#' lapply(out_3st_count_sim$emiss_naccept, function(e) e/out_3st_count_sim$input$J)
#'
#' # Finally, take the average acceptance rate by dependent variable and state:
#' lapply(out_3st_count_sim$emiss_naccept, function(e) colMeans(e/out_3st_count_sim$input$J))
#'
#'
#' @export
#'

pd_RW_emiss_count <- function(gen, emiss_scalar, emiss_w){
  if(sum(objects(gen) %in% "m") != 1 | sum(objects(gen) %in% "n_dep") != 1){
    stop("The input argument gen should contain the elements m, and n_dep.")
  }
  m <- gen$m
  n_dep <- gen$n_dep
  if(sum(sapply(emiss_scalar, is.double)) != n_dep | sum(sapply(emiss_scalar, is.matrix)) > 0 | sum(1 != sapply(emiss_scalar, length)) > 0){
    stop(paste("The elements in emiss_scalar should be a numeric vector with length 1."))
  }
  emiss_scalar <- emiss_scalar
  if(!is.list(emiss_w) | length(emiss_w) != n_dep){
    stop(paste("emiss_w should be a list containing n_dep, here", n_dep, ",elements"))
  }
  if(sum(sapply(emiss_w, is.double)) != n_dep | sum(sapply(emiss_w, is.matrix)) > 0 | sum(1 != sapply(emiss_w, length)) > 0){
    stop(paste("The elements in emiss_w should be a numeric vector with length 1."))
  }
  emiss_w <- emiss_w
  out <- list(gen = gen, emiss_scalar = emiss_scalar, emiss_w = emiss_w)
  class(out) <- append(class(out), "mHMM_pdRW_emiss")
  return(out)
}
