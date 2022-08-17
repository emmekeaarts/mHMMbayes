#' Specifying informative hyper-prior on the transition probability matrix gamma of the multilevel hidden Markov model
#'
#' \code{prior_gamma} provides a framework to manually specify an informative
#' hyper-prior on the transition probability matrix gamma, and creates an object
#' of class \code{mHMM_prior_gamma} used by the function \code{mHMM}. Note that the
#' hyper-prior distribution on the transition probabilities are on the
#' intercepts (and, if subject level covariates are used, regression
#' coefficients) of the Multinomial logit model used to accommodate the
#' multilevel framework of the data, instead of on the probabilities directly.
#' The set of hyper-prior distributions consists of a multivariate Normal
#' hyper-prior distribution on the vector of means (i.e., intercepts and
#' regression coefficients), and an Inverse Wishart hyper-prior distribution on
#' the covariance matrix.
#'
#' Estimation of the mHMM proceeds within a Bayesian context, hence a
#' hyper-prior distribution has to be defined for the group level parameters.
#' Default, non-informative priors are used unless specified otherwise by the
#' user. Each row of the transition probability matrix has its own set of
#' Multinomial logit intercepts, which are assumed to follow a multivariate
#' normal distribution. Hence, the hyper-prior distributions for the intercepts
#' consists of a multivariate Normal hyper-prior distribution on the vector of
#' means, and an Inverse Wishart hyper-prior distribution on the covariance
#' matrix. Note that only the number of states \code{m} and values of the
#' hypothesized hyper-prior mean values of the Multinomial logit intercepts have
#' to be specified by the user, default values are available for all other
#' hyper-prior distribution parameters.
#'
#' Given that the hyper-priors are specified on the intercepts of the Multinomial
#' logit model intercepts instead of on the probabilities of the transition
#' probability matrix gamma directly, specifying a hyper-prior can seem rather
#' daunting. However, see the function \code{\link{prob_to_int}} and
#' \code{\link{int_to_prob}} for translating probabilities to a set of
#' Multinomial logit intercepts and vice versa.
#'
#' Note that \code{gamma_K0}, \code{gamma_nu} and \code{gamma_V} are assumed
#' equal over the states. When the hyper-prior values for \code{gamma_K0},
#' \code{gamma_nu} and \code{gamma_V} are not manually specified, the default
#' values are as follows. \code{gamma_K0} set to 1, \code{gamma_nu} set to 3 + m
#' - 1, and the diagonal of \code{gamma_V} (i.e., the variance) set to 3 + m - 1
#' and the off-diagonal elements (i.e., the covariance) set to 0. In addition,
#' when no manual values for the hyper-prior on gamma are specified at all (that
#' is, the function \code{prior_gamma} is not used), all elements of the
#' matrices contained in \code{gamma_mu0} are set to 0 in the function
#' \code{mHMM}.
#'
#' Note that in case covariates are specified, the hyper-prior parameter values
#' of the inverse Wishart distribution on the covariance matrix remain
#' unchanged, as the estimates of the regression coefficients for the covariates
#' are fixed over subjects.
#'
#' @param m Numeric vector with length 1 denoting the number of hidden states.
#' @param n_xx_gamma Optional numeric vector with length 1 denoting the number of
#'   (level 2) covariates used to predict the transition probability matrix
#'   gamma. When omitted, the model assumes no covariates are used to predict
#'   gamma.
#' @param gamma_mu0 A list containing m matrices; one matrix for each row of the
#'   transition probability matrix gamma. Each matrix contains the hypothesized
#'   hyper-prior mean values of the intercepts of the Multinomial logit
#'   model on the transition probabilities gamma. Hence, each matrix
#'   consists of one row (when not including covariates in the model) and
#'   \code{m} - 1 columns. If covariates are used, the number of rows in each
#'   matrix in the list is equal to 1 + n_xx_gamma (i.e., the first row corresponds to
#'   the hyper-prior mean values of the intercepts, the subsequent rows
#'   correspond to the hyper-prior mean values of the regression coefficients
#'   connected to each of the covariates).
#' @param gamma_K0 Optional numeric vector with length 1 (when no covariates are
#'   used) denoting the number of hypothetical prior subjects on which the set
#'   of hyper-prior mean intercepts specified in \code{gamma_mu0} are based.
#'   When covariates are used: Numeric vector with length 1 + n_xx_gamma denoting the
#'   number of hypothetical prior subjects on which the set of intercepts (first
#'   value) and set of regression coefficients (subsequent values) are based.
#' @param gamma_nu Optional numeric vector with length 1 denoting the degrees of freedom
#'   of the hyper-prior Inverse Wishart distribution on the covariance of
#'   the Multinomial logit intercepts.
#' @param gamma_V Optional matrix of \code{m} - 1 by \code{m} - 1 containing the
#'   variance-covariance matrix of the hyper-prior Inverse Wishart distribution
#'   on the covariance of the Multinomial logit intercepts.
#'
#' @return \code{prior_gamma} returns an object of class \code{mHMM_prior_gamma},
#'   containing informative hyper-prior values for the transition probability
#'   matrix gamma of the multilevel hidden Markov model. The object is
#'   specifically created and formatted for use by the function \code{mHMM},
#'   and thoroughly checked for correct input dimensions.
#'   The object contains the following components:
#'   \describe{
#'   \item{\code{m}}{Numeric vector denoting the number of hidden states, used
#'   for checking equivalent general model properties specified under
#'   \code{prior_gamma} and \code{mHMM}.}
#'   \item{\code{gamma_mu0}}{A list containing the hypothesized hyper-prior mean
#'   values of the intercepts of the Multinomial logit model on the transition
#'   probability matrix gamma.}
#'   \item{\code{gamma_K0}}{A numeric vector denoting the number of hypothetical
#'   prior subjects on which the set of hyper-prior mean intercepts specified in
#'   \code{gamma_mu0} are based.}
#'   \item{\code{gamma_nu}}{A numeric vector denoting the degrees of freedom
#'   of the hyper-prior Inverse Wishart distribution on the covariance of the
#'   Multinomial logit intercepts.}
#'   \item{\code{gamma_V}}{A matrix containing the variance-covariance of the
#'   hyper-prior Inverse Wishart distribution on the covariance of the
#'   Multinomial logit intercepts.}
#'   \item{\code{n_xx_gamma}}{A numeric vector denoting the number of (level 2)
#'   covariates used to predict the transition probability matrix gamma. When no
#'   covariates are used, \code{n_xx_gamma} equals \code{NULL}.}
#'   }
#'
#' @seealso \code{\link{prior_emiss_cat}} for manually specifying an informative
#'   hyper-prior on the categorical emission distribution(s),
#'   \code{\link{prob_to_int}} for transforming a set of probabilities to a set
#'   of Multinomial logit regression intercepts, and \code{\link{mHMM}} for
#'   fitting a multilevel hidden Markov model.
#'
#' @examples
#' ###### Example using package example data, see ?nonverbal
#' # specifying general model properties:
#' m <- 3
#' # representing a prior belief that switching to state 3 does not occur often and
#' # state 3 has a relative short duration
#' prior_prob_gamma <- matrix(c(0.70, 0.25, 0.05,
#'                              0.25, 0.70, 0.05,
#'                              0.30, 0.30, 0.40), nrow = m, ncol = m, byrow = TRUE)
#'
#' # using the function prob_to_int to obtain intercept values for the above specified
#' # transition probability matrix gamma
#' prior_int_gamma <- prob_to_int(prior_prob_gamma)
#' gamma_mu0 <- list(matrix(prior_int_gamma[1,], nrow = 1, ncol = m-1),
#'                   matrix(prior_int_gamma[2,], nrow = 1, ncol = m-1),
#'                   matrix(prior_int_gamma[3,], nrow = 1, ncol = m-1))
#'
#' gamma_K0 <- 1
#' gamma_nu <- 5
#' gamma_V <- diag(5, m - 1)
#'
#' manual_prior_gamma <- prior_gamma(m = m, gamma_mu0 = gamma_mu0,
#'                                   gamma_K0 = gamma_K0, gamma_nu = gamma_nu,
#'                                   gamma_V = gamma_V)
#'
#'
#' # using the informative hyper-prior in a model
#' n_dep <- 4
#' q_emiss <- c(3, 2, 3, 2)
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
#' out_3st_infgamma <- mHMM(s_data = nonverbal,
#'                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                     start_val = c(list(start_TM), start_EM),
#'                     gamma_hyp_prior = manual_prior_gamma,
#'                     mcmc = list(J = 11, burn_in = 5))
#' }
#' \dontshow{
#' out_3st_infgamma <- mHMM(s_data = nonverbal,
#'                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                     start_val = c(list(start_TM), start_EM),
#'                     gamma_hyp_prior = manual_prior_gamma,
#'                     mcmc = list(J = 6, burn_in = 3))
#' }
#'
#' out_3st_infgamma
#' summary(out_3st_infgamma)
#'
#' @export
#'


prior_gamma <- function(m, gamma_mu0, gamma_K0 = NULL, gamma_nu = NULL, gamma_V = NULL, n_xx_gamma = NULL){
  if(is.null(n_xx_gamma)){
    n_xx_int <- 1
  } else {
    n_xx_int <- n_xx_gamma + 1
  }
  if(!is.list(gamma_mu0)){
    stop(paste("gamma_mu0 should be a list containing", m, "matrices; one matrix for each row of the transition probability matrix gamma."))
    }
  if(length(gamma_mu0) != m | sum(sapply(gamma_mu0, is.matrix)) != m){
    stop(paste("gamma_mu0 should be a list containing", m, "matrices; one matrix for each row of the transition probability matrix gamma."))
  }
  if(sum((m-1) == sapply(gamma_mu0, dim)[2,]) != m){
    stop(paste("Each matrix in the list gamma_mu0 should consist of m-1 , here", m-1,", columns"))
  }
  if(n_xx_int == 1 & sum(1 == sapply(gamma_mu0, dim)[1,]) != m){
    stop("According to the input paramter n_xx_gamma no covariates are used for gamma. Hence, the number of rows in each matrix of gamma_mu0 should equal 1.")
  }
  if(n_xx_int > 1 & sum(n_xx_int == sapply(gamma_mu0, dim)[1,]) != m ){
    stop(paste("According to the input paramter n_xx_gamma", n_xx_gamma, "covariates are used for gamma. Hence, the number of rows in each matrix of gamma_mu0 should equal 1 + n_xx_gamma =", 1 + n_xx_gamma))
  }
  gamma_mu0 <- gamma_mu0
  if(is.null(gamma_K0)){
    gamma_K0			<- diag(1, n_xx_int)
  } else {
    if(!is.double(gamma_K0) | length(gamma_K0) != n_xx_int){
      stop(paste("gamma_K0 should be a numeric vector with length", n_xx_int))
    }
    gamma_K0 <- diag(gamma_K0, n_xx_int)
  }
  if(is.null(gamma_nu)){
    gamma_nu <- 3 + m - 1
  } else {
    if(!is.double(gamma_nu) | length(gamma_nu) != 1){
      stop("gamma_nu should be a numeric vector with length 1")
    }
    gamma_nu <- gamma_nu
  }
  if(is.null(gamma_V)){
    gamma_V			  <- (3 + m - 1) * diag(m - 1)
  } else {
    if(!is.matrix(gamma_V)){
      stop(paste("gamma_V should be an", m-1, "by", m-1, "matrix"))
    }
    if(dim(gamma_V)[1] != (m-1) | dim(gamma_V)[2] != (m-1)){
      stop(paste("gamma_V should be an", m-1, "by", m-1, "matrix"))
    }
    gamma_V <- gamma_V
  }
  out <- list(m = m, n_xx_gamma = n_xx_gamma, gamma_mu0 = gamma_mu0, gamma_K0 = gamma_K0, gamma_nu = gamma_nu, gamma_V = gamma_V)
  class(out) <- append(class(out), "mHMM_prior_gamma")
  return(out)
}
