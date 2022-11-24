#' Proposal distribution settings RW Metropolis sampler for mHMM transition probability matrix gamma
#'
#' \code{pd_RW_gamma} provides a framework to manually specify the settings of
#' the proposal distribution of the random walk (RW) Metropolis sampler of the
#' transition probability matrix gamma of the multilevel hidden Markov model,
#' and creates on object of the class \code{mHMM_pdRW_gamma}. The RW metropolis
#' sampler is used for sampling the subject level parameter estimates relating
#' to the transition probability matrix gamma, that is, the Multinomial logistic
#' regression intercepts.
#'
#' When no manual values for the settings of the proposal distribution of the
#' random walk (RW) Metropolis sampler are specified at all (that is, the
#' function \code{pd_RW_gamma} is not used), all elements in
#' \code{gamma_int_mle0} set to 0, \code{gamma_scalar} set to 2.93 /
#' sqrt(\code{m} - 1), and \code{gamma_w} set to 0.1. See the section
#' \emph{Scaling the proposal distribution of the RW Metropolis sampler} in
#' \code{vignette("estimation-mhmm")} for details.
#'
#' Within the function \code{mHMM}, the acceptance rate of the RW metropolis
#' sampler relating to the transition probability matrix gamma can be tracked
#' using the output parameter \code{gamma_naccept}. An acceptance rate of about
#' 23\% is considered optimal when many parameters are being updated at once
#' (Gelman, Carlin, Stern & Rubin, 2014).
#'
#' @param m Numeric vector with length 1 denoting the number of hidden states.
#' @param gamma_int_mle0 A matrix with \code{m} rows and \code{m} - 1 columns
#'   denoting the starting values for the maximum likelihood (ML) estimates of
#'   the Multinomial logit regression intercepts of the transition probability
#'   matrix gamma. ML parameters to be estimated are based on the pooled data
#'   (data over all subjects).
#' @param gamma_scalar A numeric vector with length 1 denoting the scale factor
#'   \code{s}. That is, the scale of the proposal distribution is composed of a
#'   covariance matrix Sigma, which is then tuned by multiplying it by a scaling
#'   factor \code{s}^2.
#' @param gamma_w A numeric vector with length 1 denoting the weight for the
#'   overall log likelihood (i.e., log likelihood based on the pooled data over
#'   all subjects) in the fractional likelihood.
#'
#' @return \code{pd_RW_gamma} returns an object of class
#'   \code{mHMM_pdRW_gamma}, containing settings of the proposal distribution of
#'   the random walk (RW) Metropolis sampler on the transition probability
#'   matrix gamma of the multilevel hidden Markov model. The object is
#'   specifically created and formatted for use by the function \code{mHMM}, and
#'   checked for correct input dimensions. The object contains the following
#'   components:
#'    \describe{
#'   \item{\code{m}}{Numeric vector denoting the number of hidden states, used
#'   for checking equivalent general model properties specified under
#'   \code{pd_RW_gamma} and \code{mHMM}.}
#'   \item{\code{gamma_int_mle0}}{A matrix containing the starting values for
#'   the maximum likelihood (ML) estimates of the Multinomial logit regression
#'   intercepts of the transition probability matrix gamma.}
#'   \item{\code{gamma_scalar}}{A numeric vector with length 1 denoting
#'   the scale factor \code{s} of the proposal distribution.}
#'  \item{\code{gamma_w}}{A numeric vector with length 1  denoting
#'  denoting the weight for the overall log likelihood in the fractional
#'  likelihood.}
#'   }
#'
#' @references
#' \insertRef{gelman2014}{mHMMbayes}
#'
#' \insertRef{rossi2012}{mHMMbayes}
#'
#'
#' @examples
#' ###### Example using package example data, see ?nonverbal
#' # specifying general model properties:
#' m <- 3
#'
#' # specifying manual values for RW metropolis sampler on gamma
#' gamma_int_mle0 <- matrix(c( -2, -2,
#'                              2,  0,
#'                              0,  3), byrow = TRUE, nrow = m, ncol = m - 1)
#' gamma_scalar <- c(2)
#' gamma_w <- c(0.2)
#' manual_gamma_sampler <- pd_RW_gamma(m = m, gamma_int_mle0 = gamma_int_mle0,
#'                                     gamma_scalar = gamma_scalar,
#'                                     gamma_w = gamma_w)
#'
#' # specifying starting values
#' n_dep <- 4
#' q_emiss <- c(3, 2, 3, 2)
#'
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
#' out_3st_RWgamma <- mHMM(s_data = nonverbal,
#'                          gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                          start_val = c(list(start_TM), start_EM),
#'                          gamma_sampler = manual_gamma_sampler,
#'                          mcmc = list(J = 11, burn_in = 5))
#' }
#' \dontshow{
#' out_3st_RWgamma <- mHMM(s_data = nonverbal,
#'                          gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                          start_val = c(list(start_TM), start_EM),
#'                          gamma_sampler = manual_gamma_sampler,
#'                          mcmc = list(J = 6, burn_in = 3))
#' }
#'
#' out_3st_RWgamma
#' summary(out_3st_RWgamma)
#'
#' # checking acceptance rate (for illustrative purposes, in the example,
#' # J is too low for getting a fair indication)
#' J_it <- 11 - 1 # accept/reject starts at iteration 2 of MCMC algorithm
#' out_3st_RWgamma$gamma_naccept / J_it
#' # average acceptance rate over all subjects per parameter
#' apply(out_3st_RWgamma$gamma_naccept / J_it, 2, mean)
#'
#' @export
#'

pd_RW_gamma <- function(m, gamma_int_mle0, gamma_scalar, gamma_w){
  if(!is.matrix(gamma_int_mle0) | dim(gamma_int_mle0)[1] != m | dim(gamma_int_mle0)[2] != m - 1){
    stop(paste("gamma_int_mle0 should be a matrix with m - 1 columns and m rows, here", paste(paste(m, "by", m - 1), collapse = ", ")))
  }
  gamma_int_mle0 <- gamma_int_mle0
  if(!is.double(gamma_scalar) | is.matrix(gamma_scalar) | length(gamma_scalar) != 1){
    stop(paste("gamma_scalar should be a numeric vector with length 1."))
  }
  gamma_scalar <- gamma_scalar
  if(!is.double(gamma_w) | is.matrix(gamma_w) | length(gamma_w) != 1){
    stop(paste("gamma_w should be a numeric vector with length 1."))
  }
  gamma_w <- gamma_w
  out <- list(m = m, gamma_int_mle0 = gamma_int_mle0, gamma_scalar = gamma_scalar, gamma_w = gamma_w)
  class(out) <- append(class(out), "mHMM_pdRW_gamma")
  return(out)
}
