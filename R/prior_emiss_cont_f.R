#' Specifying informative hyper-prior on the continuous emission distribution(s)  of the multilevel hidden Markov model
#'
#' \code{prior_emiss_cont_f} provides a framework to manually specify an
#' informative hyper-prior on the Normal (i.e., Gaussian) emission
#' distributions. \code{prior_emiss_cont_f} creates an object of class
#' \code{mHMM_prior_emiss} used by the function \code{mHMM_f}, and additionally
#' attaches the class \code{cont} to signal use for continuous observations. The
#' set of hyper-prior distributions consists of a Normal-Inverse-Gamma
#' distribution (i.e., assuming both unknown population mean and variance
#' between subject level means) on the vector of means (i.e., intercepts and
#' regression coefficients), and an Inverse gamma distribution (i.e., assuming a
#' known mean) on each of the (fixed over subjects) emission variances.
#'
#' Estimation of the mHMM proceeds within a Bayesian context, hence a
#' hyper-prior distribution has to be defined for the group level parameters. To
#' avoid problems with 'label switching' when dealing with continuous emission
#' distribution(s) (i.e., switching of the labels of the hidden states while
#' sampling from the MCMC), the user is forced to specify hyper-prior parameter
#' values when using continuous emission distributions (i.e., default,
#' non-informative priors are not available for continuous emission
#' distributions).
#'
#' Note that \code{emiss_K0} and \code{emiss_nu} are assumed equal over the
#' states. Also note that in case covariates are specified, the hyper-prior
#' parameter values of the inverse Wishart distribution on the covariance matrix
#' remain unchanged, as the estimates of the regression coefficients for the
#' covariates are fixed over subjects.
#'
#' @inheritParams mHMM_f
#' @param emiss_mu0 A list containing \code{n_dep} matrices, i.e., one list for
#'   each dependent variable \code{k}. Each matrix contains the hypothesized
#'   hyper-prior means of the Normal emission distributions in each of the
#'   states. Hence, each matrix consists of one row (when not including
#'   covariates in the model) and \code{m} columns. If covariates are used, the
#'   number of rows in each matrix in the list is equal to 1 + n_xx (i.e., the
#'   first row corresponds to the hyper-prior means, the subsequent rows
#'   correspond to the hyper-prior values of the regression coefficients
#'   connected to each of the covariates).
#' @param emiss_K0 A list containing \code{n_dep} elements corresponding
#'   to each dependent variable \code{k}. Each element \code{k} is a
#'   numeric vector with length 1 (when no covariates are used) denoting the
#'   number of hypothetical prior subjects on which the set of hyper-prior means
#'   specified in \code{emiss_mu0} are based. When covariates are
#'   used: each element is a numeric vector with length 1 + n_xx denoting the
#'   number of hypothetical prior subjects on which the set of means (first
#'   value) and set of regression coefficients (subsequent values) are based.
#' @param emiss_V A list containing \code{n_dep} elements corresponding to each
#'   of the dependent variables \code{k}, where each element \code{k} is a
#'   vector with length \code{m} containing the hypothesized (fixed over
#'   subjects) variance of the Normal emission distributions, which are assumed
#'   to follow a Inverse Gamma hyper-prior distribution (note: here, the Inverse
#'   Gamma hyper-prior distribution is parametrized as a scaled inverse
#'   chi-squared distribution).
#' @param emiss_nu A list containing \code{n_dep} elements corresponding to each
#'   dependent variable \code{k}. Each element \code{k} is a numeric vector with
#'   length 1 denoting the degrees of freedom of the Inverse Gamma hyper-prior
#'   distribution on the (fixed over subjects) variance of the Normal emission
#'   distributions (note: here, the Inverse Gamma hyper-prior distribution is
#'   parametrized as a scaled inverse chi-squared distribution).
#'
#' @return \code{prior_emiss_cont_f} returns an object of class \code{mHMM_prior_emiss},
#'   containing informative hyper-prior values for the continuous emission
#'   distribution(s) of the multilevel hidden Markov model. The object is
#'   specifically created and formatted for use by the function \code{mHMM_f},
#'   and thoroughly checked for correct input dimensions.
#'   The object contains the following components:
#'   \describe{
#'   \item{\code{gen}}{A list containing the elements \code{m}, and \code{n_dep},
#'   used for checking equivalent general model properties
#'   specified under \code{prior_emiss_cont_f} and \code{mHMM_f}.}
#'   \item{\code{emiss_mu0}}{A lists containing the hypothesized
#'   hyper-prior means of the Normal distribution on
#'   the continuous emission probabilities.}
#'   \item{\code{emiss_K0}}{A list containing \code{n_dep} elements denoting the
#'   number of hypothetical prior subjects on which the set of hyper-prior means
#'   specified in \code{emiss_mu0} are based.}
#'   \item{\code{emiss_V}}{A list containing \code{n_dep} elements containing
#'   the variance of the Inverse Gamma hyper-prior distribution on the fixed
#'   over subjects variance of the Normal emission distribution. }
#'   \item{\code{emiss_nu}}{A list containing \code{n_dep} elements denoting the
#'   degrees of freedom of the Inverse Gamma hyper-prior distribution on the
#'   fixed over subjects variance of the Normal emission distribution. }
#'   }
#'
#' @seealso \code{\link{prior_gamma}} for manually specifying an informative
#'  hyper-prior on the transition probability matrix gamma, and
#'  \code{\link{mHMM_f}} for fitting a multilevel hidden Markov model.
#'
#' @examples
#' ###### Example using simulated data
#' # specifying general model properties:
#' m <- 3
#' n_dep <- 2
#'
#' # hypothesized hyper-prior values for the continuous emission distribution
#' manual_prior_emiss <- prior_emiss_cont_f(
#'                         gen = list(m = m, n_dep = n_dep),
#'                         emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
#'                                          matrix(c(7, 8, 18), nrow = 1)),
#'                         emiss_K0 = list(1, 1),
#'                         emiss_V =  list(rep(5^2, m), rep(0.5^2, m)),
#'                         emiss_nu = list(1, 1))
#'
#' # to use the informative priors in a model, simulate multivariate continuous data
#' n_t     <- 100
#' n       <- 10
#'
#' gamma   <- matrix(c(0.8, 0.1, 0.1,
#'                     0.2, 0.7, 0.1,
#'                     0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
#'
#' emiss_distr <- list(matrix(c( 50, 10,
#'                               100, 10,
#'                               150, 10), nrow = m, byrow = TRUE),
#'                     matrix(c(5, 2,
#'                              10, 5,
#'                              20, 3), nrow = m, byrow = TRUE))
#'
#' data_cont <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous',
#'                       gen = list(m = m, n_dep = n_dep),
#'                       gamma = gamma, emiss_distr = emiss_distr,
#'                       var_gamma = .1, var_emiss = c(5^2, 0.2^2))
#'
#' # using the informative hyper-prior in a model
#' # Note that for reasons of running time, J is set at a ridiculous low value.
#' # One would typically use a number of iterations J of at least 1000,
#' # and a burn_in of 200.
#' out_3st_cont_sim_infemiss <- mHMM_f(s_data = data_cont$obs,
#'                     data_distr = "continuous",
#'                     gen = list(m = m, n_dep = n_dep),
#'                     start_val = c(list(gamma), emiss_distr),
#'                     emiss_hyp_prior = manual_prior_emiss,
#'                     mcmc = list(J = 11, burn_in = 5))
#'
#' out_3st_cont_sim_infemiss
#' summary(out_3st_cont_sim_infemiss)
#'
#' @export
#'


prior_emiss_cont_f <- function(gen, emiss_mu0, emiss_K0, emiss_V, emiss_nu){
  if(sum(objects(gen) %in% "m") != 1 | sum(objects(gen) %in% "n_dep") != 1){
    stop("The input argument gen should contain the elements m and n_dep")
  }
  m <- gen$m
  n_dep <- gen$n_dep

  #### checking emiss_mu0 ####
  if(!is.list(emiss_mu0)){
    stop(paste("emiss_mu0 should be a list containing", n_dep, "matrices; one matrix for each dependent variable."))
  }
  if(length(emiss_mu0) != n_dep ){
    stop(paste("emiss_mu0 should be a list containing", n_dep, "matrices; one matrix for each dependent variable."))
  }
  if(any(!sapply(emiss_mu0, is.matrix))){
    stop(paste("emiss_mu0 should be a list containing", n_dep, "matrices; one matrix for each dependent variable."))
  }
  if(sum(m == sapply(emiss_mu0, dim)[2,]) != n_dep){
    stop(paste("The matrix relating to dependent variable", k, "of the input argument emiss_mu0 should consist of m, here", m, ", columns."))
  }

  #### checking emiss_K0 ####
  if(!is.list(emiss_K0) | length(emiss_K0) != n_dep){
    stop(paste("emiss_K0 should be a list containing n_dep, here", n_dep,", elements. "))
  }
  if(sum(sapply(emiss_K0, is.double)) != n_dep | sum(sapply(emiss_K0, is.matrix)) > 0 ){
    stop(paste("Each of the n_dep elements within the list emiss_K0 should be a numeric vector with lengths", paste(1, collapse = ", "), "."))
  }
  emiss_K0_list     <- rep(list(NULL), n_dep)
  for(k in 1:n_dep){
    emiss_K0_list[[k]]			<- diag(emiss_K0[[k]], 1)
  }
  emiss_K0 <- emiss_K0_list

  #### checking emiss_nu ####
  if(!is.list(emiss_nu)){
     stop(paste("emiss_nu should be a list containing n_dep, here", n_dep,", elements, where each element is a numeric vector with length 1 "))
  }
  if(length(emiss_nu) != n_dep | sum(rep(1, n_dep) == sapply(emiss_nu, length)) != n_dep |
     sum(sapply(emiss_nu, is.double)) != n_dep | sum(sapply(emiss_nu, is.matrix)) > 0){
    stop(paste("emiss_nu should be a list containing n_dep, here", n_dep,", elements, where each element is a numeric vector with length 1 "))
  }

  #### checking emiss_V ####
  if(!is.list(emiss_V)){
     stop(paste("emiss_V should be a list containing n_dep, here", n_dep,", elements, where each element is vector with lenght m, here", m, "."))
  }
  if(length(emiss_V) != n_dep | sum(m == sapply(emiss_V, length)) != n_dep){
    stop(paste("emiss_V should be a list containing n_dep, here", n_dep,", elements, where each element is vector with lenght m, here", m, "."))
  }



  #### return output ####
  out <- list(gen = gen, emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0, emiss_nu = emiss_nu, emiss_V = emiss_V)
  class(out) <- append(class(out), c("mHMM_prior_emiss", "cont_f"))
  return(out)
}
