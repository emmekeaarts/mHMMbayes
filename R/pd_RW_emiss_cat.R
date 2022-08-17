#' Proposal distribution settings RW Metropolis sampler for mHMM categorical emission distribution(s)
#'
#' \code{pd_RW_emiss_cat} provides a framework to manually specify the
#' settings of the proposal distribution of the random walk (RW) Metropolis
#' sampler of the emission distribution(s) of the multilevel hidden Markov
#' model, and creates on object of the class \code{mHMM_pdRW_emiss}. The RW
#' metropolis sampler is used for sampling the subject level parameter estimates
#' relating to the emission distributions of the dependent variables \code{k},
#' that is, the Multinomial logistic regression intercepts.
#'
#'
#' When no manual values for the settings of the proposal distribution of the
#' random walk (RW) Metropolis sampler are specified at all (that is, the
#' function \code{pd_RW_emiss_cat} is not used), all elements in
#' \code{emiss_int_mle0} set to 0, \code{emiss_scalar} set to 2.93 /
#' sqrt(\code{q_emiss[k]} - 1), and \code{emiss_w} set to 0.1. See the section
#' \emph{Scaling the proposal distribution of the RW Metropolis sampler} in
#' \code{vignette("estimation-mhmm")} for details.
#'
#' Within the function \code{mHMM}, the acceptance rate of the RW metropolis
#' sampler relating to the emission distribution(s) can be tracked using the
#' output parameter \code{emiss_naccept}. An acceptance rate of about 23\% is
#' considered optimal when many parameters are being updated at once (Gelman,
#' Carlin, Stern & Rubin, 2014).
#'
#' @inheritParams mHMM
#' @param emiss_int_mle0 A list containing \code{n_dep} elements corresponding
#'   to each of the dependent variables \code{k}, where each element is a matrix
#'   with \code{m} rows and \code{q_emiss[k]} - 1 columns denoting the starting
#'   values for the maximum likelihood (ML) estimates of the Multinomial logit
#'   regression intercepts of the emission distribution(s). ML parameters to be
#'   estimated are based on the pooled data (data over all subjects).
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
#' @return \code{pd_RW_emiss_cat} returns an object of class
#'   \code{mHMM_pdRW_emiss}, containing settings of the proposal distribution of
#'   the random walk (RW) Metropolis sampler on the categorical emission
#'   distribution(s) of the multilevel hidden Markov model. The object is
#'   specifically created and formatted for use by the function \code{mHMM}, and
#'   checked for correct input dimensions. The object contains the following
#'   components:
#'    \describe{
#'   \item{\code{gen}}{A list containing the elements \code{m}, \code{n_dep},
#'   and \code{q_emiss}, used for checking equivalent general model properties
#'   specified under \code{pd_RW_emiss_cat} and \code{mHMM}.}
#'   \item{\code{emiss_int_mle0}}{A list containing \code{n_dep} elements, where
#'   each element is a matrix  containing the starting values for the maximum
#'   likelihood (ML) estimates of the Multinomial logit regression intercepts of
#'   the emission distribution(s).}
#'   \item{\code{emiss_scalar}}{A list containing \code{n_dep} elements denoting
#'   the scale factor \code{s} of the proposal distribution.}
#'  \item{\code{emiss_w}}{A list containing \code{n_dep} elements denoting
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
#' n_dep <- 4
#' q_emiss <- c(3, 2, 3, 2)
#'
#' # specifying manual values for RW metropolis sampler on emission distribtutions
#' emiss_int_mle0 <- list(matrix(c( 2,  0,
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
#' emiss_scalar <- list(c(2), c(3), c(2), c(3))
#' emiss_w <- rep(list(c(0.2)), n_dep)
#' manual_emiss_sampler <- pd_RW_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                                         emiss_int_mle0 = emiss_int_mle0,
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
#' }
#' \dontshow{
#' out_3st_RWemiss <- mHMM(s_data = nonverbal,
#'                          gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                          start_val = c(list(start_TM), start_EM),
#'                          emiss_sampler = manual_emiss_sampler,
#'                          mcmc = list(J = 6, burn_in = 3))
#' }
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
#' t(sapply(RW_emiss_accept, apply, MARGIN = 2, mean, simplyfy = FALSE))
#'
#' @export
#'

pd_RW_emiss_cat <- function(gen, emiss_int_mle0, emiss_scalar, emiss_w){
  if(sum(objects(gen) %in% "m") != 1 | sum(objects(gen) %in% "n_dep") != 1 | sum(objects(gen) %in% "q_emiss") != 1){
    stop("The input argument gen should contain the elements m, n_dep and q_emiss.")
  }
  m <- gen$m
  n_dep <- gen$n_dep
  q_emiss <- gen$q_emiss
  if(length(q_emiss) != n_dep){
    stop("The lenght of q_emiss specifying the number of output categories for each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
  }
  if(!is.list(emiss_int_mle0) | length(emiss_int_mle0) != n_dep){
    stop(paste("emiss_int_mle0 should be a list containing n_dep, here", n_dep, ",elements"))
  }
  if(sum(sapply(emiss_int_mle0, is.matrix)) != n_dep | sum(q_emiss - 1 != sapply(emiss_int_mle0, dim)[2,]) > 0 | sum(m != sapply(emiss_int_mle0, dim)[1,]) > 0){
    stop(paste("The elements in emiss_int_mle0 should be a matrix with q_emiss - 1 columns and m rows, here", paste(paste(m, "by", q_emiss - 1), collapse = ", ")))
  }
  emiss_int_mle0 <- emiss_int_mle0
  if(!is.list(emiss_scalar) | length(emiss_scalar) != n_dep){
    stop(paste("emiss_scalar should be a list containing n_dep, here", n_dep, ",elements"))
  }
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
  out <- list(gen = gen, emiss_int_mle0 = emiss_int_mle0, emiss_scalar = emiss_scalar, emiss_w = emiss_w)
  class(out) <- append(class(out), "mHMM_pdRW_emiss")
  return(out)
}
