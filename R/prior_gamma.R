#' Specifying informative hyper-prior on the transition probability matrix gamma of the multilevel hidden Markov model
#'
#' \code(prior_gamma) provides a framework to manually specify an informative
#' hyper-prior on the transition probability matrix gamma, and creates an object
#' of class mHMM_prior_gamma used by the function \code{mHMM}. Estimation of the
#' mHMM proceeds within a Bayesian context, hence a hyper-prior distribution has
#' to be defined for the group level parameters. Default, non-informative priors
#' are used unless specified otherwise by the user. Note that the hyper-prior
#' distribution on the transition probabilities are on the intercepts (and, if
#' subject level covariates are used, regression coefficients) of the
#' Multinomial regression model used to accommodate the multilevel framework of
#' the data, instead of on the probabilities directly. See the function
#' \code{\link{prob_to_int}} and \code{\link{int_to_prob}} for translating
#' probabilities to a set of Multinomial regression intercepts and vice versa.
#' Each row of the transition probability matrix has its own set of Multinomial
#' regression intercepts, which are assumed to follow a multivariate normal
#' distribution. Hence, the hyper-prior distributions for the intercepts
#' consists of a multivariate Normal hyper-prior distribution on the vector of
#' means, and an Inverse Wishart hyper-prior distribution on the covariance
#' matrix.
#'
#' Note that \code{gamma_K0}, \code{gamma_nu} and \code{gamma_V} are assumed
#' equal over the states. When the hyper-prior on the transition probability
#' matrix gamma is not manually specified, the default values are as follows.
#' All elements of the matrices contained in \code{gamma_mu0} are set to 0,
#' \code{gamma_K0} set to 1, \code{gamma_nu} set to 3 + m - 1, and the diagonal
#' of \code{gamma_V} (i.e., the variance) set to 3 + m - 1 and the off-diagonal
#' elements (i.e., the covariance) set to 0.
#'
#' Note that in case covariates are specified, the hyper-prior parameter values
#' of the inverse Wishart distribution on the covariance matrix remain
#' unchanged, as the estimates of the regression coefficients for the covariates
#' are fixed over subjects.
#'
#' @param m Numeric vector with length 1 denoting the number of hidden states.
#' @param n_xx Optional numeric vector with length 1 denoting the number of
#'   (level 2) covariates used to predict the transition probability matrix
#'   gamma. When omitted, the model assumes no covariates are used.
#' @param gamma_mu0 A list containing m matrices; one matrix for each row of the
#'   transition probability matrix gamma. Each matrix contains the hypothesized
#'   mean values of the intercepts. Hence, each matrix consists of one row (when
#'   not including covariates in the model) and \code{m} - 1 columns. If
#'   covariates are used, the number of rows in each matrix in the list is
#'   equal to 1 + n_xx (i.e., the first row corresponds to the hyper-prior mean
#'   values of the intercepts, the subsequent rows correspond to the hyper-prior
#'   mean values of the regression coefficients connected to each of the
#'   covariates).
#' @param gamma_K0 Optional numeric vector with length 1 (when no covariates are
#'   used) denoting the number of hypothetical prior subjects on which the set
#'   of mean intercepts specified in \code{gamma_mu0} are based. When covariates
#'   are used: Numeric vector with length 1 + n_xx denoting the number of
#'   hypothetical prior subjects on which the set of intercepts (first value)
#'   and set of regression coefficients (subsequent values) are based. When
#'   omitted, the default value of 1 is used for all elements of the numeric
#'   vector.
#' @param gamma_nu Numeric vector with length 1 denoting the degrees of
#'  freedom of the Inverse Wishart distribution
#' @param gamma_V Matrix of \code{m} - 1 by \code{m} - 1 containing the
#'  hypothesized variance-covariance matrix between the set of intercepts.
#'
#' @examples
#' # specifying general model properties:
#' m <- 3
#' # representing a prior belief that switching to state 3 does not occur often and state 3 has a relative short duration
#' prior_prob_gamma <- matrix(c(0.70, 0.25, 0.05,
#'                              0.25, 0.70, 0.05,
#'                              0.30, 0.30, 0.40), nrow = m, ncol = m, byrow = TRUE)
#'
#' # using the function prob_to_int to obtain intercept values for the above specified transition probability matrix gamma
#' prior_int_gamma <- prob_to_int(prior_prob_gamma)
#' gamma_mu0 <- list(matrix(prior_int_gamma[1,], nrow = 1, ncol = m-1),
#'                   matrix(prior_int_gamma[2,], nrow = 1, ncol = m-1),
#'                   matrix(prior_int_gamma[3,], nrow = 1, ncol = m-1))
#'
#' gamma_K0 <- 1
#' gamma_nu <- 5
#' gamma_V <- diag(5, m)
#'
#' manual_prior_gamma <- prior_gamma(m = m, gamma_mu0 = gamma_mu0,
#'                                   gamma_K0 = gamma_K0, gamma_nu = gamma_nu, gamma_V = gamma_V)
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
#' out_3st_infgamma <- mHMM(s_data = nonverbal,
#'                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                     start_val = c(list(start_TM), start_EM),
#'                     gamma_hyp_prior = manual_prior_gamma,
#'                     mcmc = list(J = 11, burn_in = 5))
#'
#' out_3st_infgamma
#' summary(out_3st_infgamma)
#'
#' @export
#'
#' INCLUDE EXAMPLE WHERE WE USE prob_to_int directly to obtain a

prior_gamma <- function(m, n_xx = NULL, gamma_mu0, gamma_K0 = NULL, gamma_nu = NULL, gamma_V = NULL){
  if(is.null(n_xx)){
    n_xx_int <- 1
  } else {
    n_xx_int <- n_xx + 1
  }
  if(!is.list(gamma_mu0) | length(gamma_mu0) != m){

  }
  gamma_mu0 <- gamma_mu0
  if(is.null(gamma_K0)){
    gamma_K0			<- diag(1, n_xx_int)
  } else {
    gamma_K0 <- diag(gamma_K0, n_xx_int)
  }
  if(is.null(gamma_nu)){
    gamma_nu <- 3 + m - 1
  } else {
    gamma_nu <- gamma_nu
  }
  if(is.null(gamma_V)){
    gamma_V			  <- (3 + m - 1) * diag(m - 1)
  } else {
    gamma_V <- gamma_V
  }
  out <- list(m = m, n_xx = n_xx, gamma_mu0 = gamma_mu0, gamma_K0 = gamma_K0, gamma_nu = gamma_nu, gamma_V = gamma_V)
  class(out) <- append(class(out), "mHMM_prior_gamma")
  return(out)
}
