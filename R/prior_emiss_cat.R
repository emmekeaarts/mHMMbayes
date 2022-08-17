#' Specifying informative hyper-prior on the categorical emission distribution(s)  of the multilevel hidden Markov model
#'
#' \code{prior_emiss_cat} provides a framework to manually specify an informative
#' hyper-prior on the categorical emission distribution(s), and creates an object
#' of class \code{mHMM_prior_emiss} used by the function \code{mHMM}. Note that the
#' hyper-prior distribution on the categorical emission probabilities are on the
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
#' user. For each dependent variable, each row of the categorical emission
#' probability matrix (i.e., the probability to observe each category (columns)
#' within each of the states (rows)) has its own set of Multinomial logit
#' intercepts, which are assumed to follow a multivariate normal distribution.
#' Hence, the hyper-prior distributions for the intercepts consists of a
#' multivariate Normal hyper-prior distribution on the vector of means, and an
#' Inverse Wishart hyper-prior distribution on the covariance matrix. Note that
#' only the general model properties (number of states \code{m}, number of
#' dependent variables \code{n_dep}, and number of observed categories for each
#' dependent variable \code{q_emiss}) and values of the hypothesized hyper-prior
#' mean values of the Multinomial logit intercepts have to be specified by the
#' user, default values are available for all other hyper-prior distribution
#' parameters.
#'
#' Given that the hyper-priors are specified on the intercepts of the Multinomial
#' logit model intercepts instead of on the categorical emission
#' probabilities directly, specifying a hyper-prior can seem rather
#' daunting. However, see the function \code{\link{prob_to_int}} and
#' \code{\link{int_to_prob}} for translating probabilities to a set of
#' Multinomial logit intercepts and vice versa.
#'
#' Note that \code{emiss_K0}, \code{emiss_nu} and \code{emiss_V} are assumed
#' equal over the states.  When the hyper-prior values for \code{emiss_K0},
#' \code{emiss_nu} and \code{emiss_V} are not manually specified, the default
#' values are as follows. \code{emiss_K0} set to 1, \code{emiss_nu} set to 3 +
#' \code{q_emiss[k]} - 1, and the diagonal of \code{gamma_V} (i.e., the
#' variance) set to 3 + \code{q_emiss[k]} - 1 and the off-diagonal elements
#' (i.e., the covariance) set to 0. In addition, when no manual values for the
#' hyper-prior on the categorical emission distribution are specified at all
#' (that is, the function \code{prior_emiss_cat} is not used), all elements of
#' the matrices contained in \code{emiss_mu0} are set to 0 in the function
#' \code{mHMM}.
#'
#' Note that in case covariates are specified, the hyper-prior parameter values
#' of the inverse Wishart distribution on the covariance matrix remain
#' unchanged, as the estimates of the regression coefficients for the covariates
#' are fixed over subjects.
#'
#' @inheritParams mHMM
#' @param n_xx_emiss Optional numeric vector with length \code{n_dep} denoting
#'   the number of (level 2) covariates used to predict the emission
#'   distribution of each of the dependent variables \code{k}. When omitted, the
#'   model assumes no covariates are used to predict the emission
#'   distribution(s).
#' @param emiss_mu0 A list of lists: \code{emiss_mu0} contains \code{n_dep}
#'   lists, i.e., one list for each dependent variable \code{k}. Each list
#'   \code{k} contains \code{m} matrices; one matrix for each set of emission
#'   probabilities within a state. The matrices contain the hypothesized
#'   hyper-prior mean values of the intercepts of the Multinomial logit model on
#'   the categorical emission probabilities. Hence, each matrix consists of one
#'   row (when not including covariates in the model) and \code{q_emiss[k]} - 1
#'   columns. If covariates are used, the number of rows in each matrix in the
#'   list is equal to 1 + n_xx (i.e., the first row corresponds to the
#'   hyper-prior mean values of the intercepts, the subsequent rows correspond
#'   to the hyper-prior mean values of the regression coefficients connected to
#'   each of the covariates).
#' @param emiss_K0 Optional list containing \code{n_dep} elements corresponding
#'   to each dependent variable \code{k}. Each element \code{k} is a
#'   numeric vector with length 1 (when no covariates are used) denoting the
#'   number of hypothetical prior subjects on which the set of hyper-prior mean
#'   intercepts specified in \code{emiss_mu0} are based. When covariates are
#'   used: each element is a numeric vector with length 1 + n_xx denoting the
#'   number of hypothetical prior subjects on which the set of intercepts (first
#'   value) and set of regression coefficients (subsequent values) are based.
#' @param emiss_nu Optional list containing \code{n_dep} elements corresponding
#'   to each dependent variable \code{k}. Each element \code{k} is
#'   a numeric vector with length 1 denoting the degrees of freedom of the
#'   hyper-prior Inverse Wishart distribution on the covariance of the
#'   Multinomial logit intercepts.
#' @param emiss_V Optional list containing \code{n_dep} elements corresponding
#'   to each dependent variable \code{k}, where each element \code{k} is a
#'   matrix of \code{q_emiss[k]} - 1 by \code{q_emiss[k]} - 1 containing the
#'   variance-covariance of the hyper-prior Inverse Wishart distribution on the
#'   covariance of the Multinomial logit intercepts.
#'
#' @return \code{prior_emiss_cat} returns an object of class \code{mHMM_prior_emiss},
#'   containing informative hyper-prior values for the categorical emission
#'   distribution(s) of the multilevel hidden Markov model. The object is
#'   specifically created and formatted for use by the function \code{mHMM},
#'   and thoroughly checked for correct input dimensions.
#'   The object contains the following components:
#'   \describe{
#'   \item{\code{gen}}{A list containing the elements \code{m}, \code{n_dep},
#'   and \code{q_emiss}, used for checking equivalent general model properties
#'   specified under \code{prior_emiss_cat} and \code{mHMM}.}
#'   \item{\code{emiss_mu0}}{A list of lists containing the hypothesized
#'   hyper-prior mean values of the intercepts of the Multinomial logit model on
#'   the categorical emission probabilities.}
#'   \item{\code{emiss_K0}}{A list containing \code{n_dep} elements denoting the
#'   number of hypothetical prior subjects on which the set of hyper-prior mean
#'   intercepts specified in \code{emiss_mu0} are based.}
#'   \item{\code{emiss_nu}}{A list containing \code{n_dep} elements denoting the
#'   degrees of freedom of the hyper-prior Inverse Wishart distribution on
#'   the covariance of the Multinomial logit intercepts.}
#'   \item{\code{emiss_V}}{A list containing \code{n_dep} elements containing
#'   the variance-covariance of the hyper-prior Inverse Wishart distribution on
#'   the covariance of the Multinomial logit intercepts.}
#'   \item{\code{n_xx_emiss}}{A numeric vector denoting the number of (level 2)
#'   covariates used to predict the emission distribution of each of the
#'   dependent variables. When no covariates are used, \code{n_xx_emiss} equals
#'   \code{NULL}.}
#'   }
#'
#' @seealso \code{\link{prior_gamma}} for manually specifying an informative
#'   hyper-prior on the transition probability matrix gamma,
#'   \code{\link{prob_to_int}} for transforming a set of probabilities to a set
#'   of Multinomial logit regression intercepts, and \code{\link{mHMM}} for
#'   fitting a multilevel hidden Markov model.
#'
#' @examples
#' ###### Example using package example data, see ?nonverbal
#' # specifying general model properties:
#' m <- 3
#' n_dep <- 4
#' q_emiss <- c(3, 2, 3, 2)
#'
#' # hypothesized mean emission probabilities
#' prior_prob_emiss_cat <- list(matrix(c(0.10, 0.80, 0.10,
#'                                       0.80, 0.10, 0.10,
#'                                       0.40, 0.40, 0.20), byrow = TRUE,
#'                                     nrow = m, ncol = q_emiss[1]), # vocalizing patient,
#'                                     # prior belief: state 1 - much talking, state 2 -
#'                                     # no talking, state 3 - mixed
#'                              matrix(c(0.30, 0.70,
#'                                       0.30, 0.70,
#'                                       0.30, 0.70), byrow = TRUE, nrow = m,
#'                                     ncol = q_emiss[2]), # looking patient
#'                                     # prior belief: all 3 states show frequent looking
#'                                     # behavior
#'                              matrix(c(0.80, 0.10, 0.10,
#'                                       0.10, 0.80, 0.10,
#'                                       0.40, 0.40, 0.20), byrow = TRUE,
#'                                     nrow = m, ncol = q_emiss[3]), # vocalizing therapist
#'                                     # prior belief: state 1 - no talking, state 2 -
#'                                     # frequent talking, state 3 - mixed
#'                              matrix(c(0.30, 0.70,
#'                                       0.30, 0.70,
#'                                       0.30, 0.70), byrow = TRUE, nrow = m,
#'                                     ncol = q_emiss[4])) # looking therapist
#'                                     # prior belief: all 3 states show frequent looking
#'                                     # behavior
#'
#' # using the function prob_to_int to obtain intercept values for the above specified
#' # categorical emission distributions
#' prior_int_emiss <- sapply(prior_prob_emiss_cat, prob_to_int)
#' emiss_mu0 <- rep(list(vector(mode = "list", length = m)), n_dep)
#' for(k in 1:n_dep){
#'   for(i in 1:m){
#'   emiss_mu0[[k]][[i]] <- matrix(prior_int_emiss[[k]][i,], nrow = 1)
#'   }
#' }
#'
#' emiss_K0 <- rep(list(c(1)), n_dep)
#' emiss_nu <- list(c(5), c(4), c(5), c(4))
#' emiss_V <- list(diag(5, q_emiss[1] - 1),
#'                 diag(4, q_emiss[2] - 1),
#'                 diag(5, q_emiss[3] - 1),
#'                 diag(4, q_emiss[4] - 1))
#'
#' manual_prior_emiss <- prior_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                                   emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
#'                                   emiss_nu = emiss_nu, emiss_V = emiss_V)
#'
#'
#' # using the informative hyper-prior in a model
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
#' out_3st_infemiss <- mHMM(s_data = nonverbal,
#'                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                     start_val = c(list(start_TM), start_EM),
#'                     emiss_hyp_prior = manual_prior_emiss,
#'                     mcmc = list(J = 11, burn_in = 5))
#' }
#' \dontshow{
#' # executable in < 5 sec together with the examples above
#' out_3st_infemiss <- mHMM(s_data = nonverbal,
#'                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                     start_val = c(list(start_TM), start_EM),
#'                     emiss_hyp_prior = manual_prior_emiss,
#'                     mcmc = list(J = 6, burn_in = 3))
#' }
#'
#' out_3st_infemiss
#' summary(out_3st_infemiss)
#'
#' @export
#'


prior_emiss_cat <- function(gen, emiss_mu0, emiss_K0 = NULL, emiss_nu = NULL, emiss_V = NULL, n_xx_emiss = NULL){
  if(sum(objects(gen) %in% "m") != 1 | sum(objects(gen) %in% "n_dep") != 1 | sum(objects(gen) %in% "q_emiss") != 1){
    stop("The input argument gen should contain the elements m, n_dep and q_emiss.")
  }
  m <- gen$m
  n_dep <- gen$n_dep
  q_emiss <- gen$q_emiss
  if(is.null(n_xx_emiss)){
    n_xx_int <- rep(1, n_dep)
  } else {
    if(length(n_xx_emiss) != n_dep){
      stop(paste("n_xx_emiss should be a numeric vector with length n_dep, here", n_dep, "."))
    }
    n_xx_int <- n_xx_emiss + 1
  }
  if(length(q_emiss) != n_dep){
    stop("The lenght of q_emiss specifying the number of output categories for each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
  }
  if(!is.list(emiss_mu0)){
    stop(paste("emiss_mu0 should be a list of lists containing", n_dep, "lists; one list for each dependent variable."))
  }
  if(length(emiss_mu0) != n_dep | sum(sapply(emiss_mu0, is.list)) != n_dep){
    stop(paste("emiss_mu0 should be a list of lists containing", n_dep, "lists; one list for each dependent variable."))
  }
  if(sum(m == sapply(emiss_mu0, length)) != n_dep){
    stop(paste("Each list k in the input argument emiss_mu0 should consist of m , here", m,", elements"))
  }
  for(k in 1:n_dep){
    if(sum(sapply(emiss_mu0[[k]], is.matrix)) != m){
      stop(paste("For the input argument emiss_mu0, each element within list k should be a matrix. This does not hold for the lists relating to dependent variable", k))
    }
    if(sum((q_emiss[k] - 1) == sapply(emiss_mu0[[k]], dim)[2,]) != m){
      stop(paste("For the list relating to dependent variable", k, "of the input argument emiss_mu0, each element should consists of a matrix with q_emiss[k] - 1, here", q_emiss[k] - 1, ", columns."))
    }
    if(n_xx_int[k] == 1 & sum(1 == sapply(emiss_mu0[[k]], dim)[1,]) != m){
      stop(paste("According to the input paramter n_xx_emiss no covariates are used to predict the emission distribution of dependent variable", k, ". Hence, within input argument emiss_mu0, in the list relating to dependent variable", k, ",each element should be a matrix with 1 row."))
    }
    if(n_xx_int[k] > 1 & sum(n_xx_int[k] == sapply(emiss_mu0[[k]], dim)[1,]) != m ){
      stop(paste("According to the input paramter n_xx_emiss", n_xx_emiss[k], "covariates are used to predict the emission distribution of dependent variable", k, ". Hence, within input argument emiss_mu0, in the list relating to dependent variable", k, ", each element should be a matrix with 1 + n_xx_emiss =", 1 + n_xx_emiss[k], "rows."))
    }
  }
  emiss_mu0 <- emiss_mu0
  if(is.null(emiss_K0)){
    emiss_K0     <- rep(list(NULL), n_dep)
    for(k in 1:n_dep){
      emiss_K0[[k]]			<- diag(1, n_xx_int[k])
    }
  } else {
    if(!is.list(emiss_K0) | length(emiss_K0) != n_dep){
      stop(paste("emiss_K0 should be a list containing n_dep, here", n_dep,", elements. "))
    }
    if(sum(sapply(emiss_K0, is.double)) != n_dep | sum(sapply(emiss_K0, is.matrix)) > 0 | sum(n_xx_int == sapply(emiss_K0, length)) != n_dep){
      stop(paste("Each of the n_dep elements within the list emiss_K0 should be a numeric vector with lengths", paste(n_xx_int, collapse = ", "), "."))
    }
    emiss_K0_list     <- rep(list(NULL), n_dep)
    for(k in 1:n_dep){
      emiss_K0_list[[k]]			<- diag(emiss_K0[[k]], n_xx_int[k])
    }
    emiss_K0 <- emiss_K0_list
  }
  if(is.null(emiss_nu)){
    emiss_nu	    <- rep(list(NULL), n_dep)
    for(k in 1:n_dep){
      emiss_nu[[k]]		<- 3 + q_emiss[k] - 1
    }
  } else {
    if(!is.list(emiss_nu)){
      stop(paste("emiss_nu should be a list containing n_dep, here", n_dep,", elements, where each element is a numeric vector with length 1 "))
    }
    if(length(emiss_nu) != n_dep | sum(rep(1, n_dep) == sapply(emiss_nu, length)) != n_dep |
       sum(sapply(emiss_nu, is.double)) != n_dep | sum(sapply(emiss_nu, is.matrix)) > 0){
      stop(paste("emiss_nu should be a list containing n_dep, here", n_dep,", elements, where each element is a numeric vector with length 1 "))
    }
    emiss_nu <- emiss_nu
  }
  if(is.null(emiss_V)){
    emiss_V	    <- rep(list(NULL), n_dep)
    for(k in 1:n_dep){
      emiss_V[[k]]		  <- (3 + q_emiss[k] - 1) * diag(q_emiss[k] - 1)
    }
  } else {
    if(!is.list(emiss_V)){
      stop(paste("emiss_V should be a list containing n_dep, here", n_dep,", elements, where each element is a matrix with the following dimensions: ", paste(paste(q_emiss - 1, "by", q_emiss - 1), collapse = ", "), "."))
    }
    if(length(emiss_V) != n_dep | sum(sapply(emiss_V, is.matrix)) != n_dep | sum(q_emiss-1 == sapply(emiss_V, dim)[1,]) != n_dep | sum(q_emiss-1 == sapply(emiss_V, dim)[2,]) != n_dep){
      stop(paste("emiss_V should be a list containing n_dep, here", n_dep,", elements, where each element is a matrix with the following dimensions: ", paste(paste(q_emiss - 1, "by", q_emiss - 1), collapse = ", "), "."))
    }
    emiss_V <- emiss_V
  }
  out <- list(gen = gen, emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0, emiss_nu = emiss_nu, emiss_V = emiss_V, n_xx_emiss = n_xx_emiss)
  class(out) <- append(class(out), "mHMM_prior_emiss")
  return(out)
}
