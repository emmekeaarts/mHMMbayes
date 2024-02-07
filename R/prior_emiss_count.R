#' Specifying informative hyper-prior on the continuous emission distribution(s)  of the multilevel hidden Markov model
#'
#' \code{prior_emiss_count} provides a framework to manually specify an
#' informative hyper-prior on the Poisson-lognormal emission
#' distributions. \code{prior_emiss_count} creates an object of class
#' \code{mHMM_prior_emiss} used by the function \code{mHMM}, and additionally
#' attaches the class \code{count} to signal use for count observations. The
#' set of hyper-prior distributions consists of a lognormal-Inverse-Gamma
#' distribution (i.e., assuming both unknown population mean and variance
#' between subject level means) on the vector of poisson means (i.e.,
#' intercepts and regression coefficients).
#'
#' Estimation of the mHMM proceeds within a Bayesian context, hence a
#' hyper-prior distribution has to be defined for the group level parameters. To
#' avoid problems with 'label switching' when dealing with continuous emission
#' distribution(s) (i.e., switching of the labels of the hidden states while
#' sampling from the MCMC), the user is forced to specify hyper-prior parameter
#' values when using count emission distributions (i.e., default,
#' non-informative priors are not available for count emission
#' distributions).
#'
#' Note that \code{emiss_K0} and \code{emiss_nu} are assumed equal over the
#' states. Also note that in case covariates are specified, the hyper-prior
#' parameter values of the inverse Wishart distribution on the covariance matrix
#' remain unchanged, as the estimates of the regression coefficients for the
#' covariates are fixed over subjects.
#'
#' Also note that for simplicity the hyper-prior means and variances of the
#' lognormal distribution, \code{emiss_mu0} and \code{emiss_V}, by default
#' have to be specified in the natural (positive real numbers) scale and not in
#' the logarithmic scale. \code{prior_emiss_count} returns the corresponding
#' values of the parameters on the logarithmic scale. If the user wants to
#' manually specify these values on the logarithmic scale, please set the
#' argument \code{log_scale} to \code{TRUE} in \code{prior_emiss_count}. If
#' covariates are used to predict the emission distribution, then the
#' logarithmic scale should be used for the inputs \code{emiss_mu0} and
#' \code{emiss_V}, and set \code{log_scale = TRUE}.
#'
#' @inheritParams mHMM
#' @param n_xx_emiss Optional numeric vector with length \code{n_dep} denoting
#'   the number of (level 2) covariates used to predict the emission
#'   distribution of each of the dependent variables \code{k}. When omitted, the
#'   model assumes no covariates are used to predict the emission
#'   distribution(s).
#' @param emiss_mu0 A list containing \code{n_dep} matrices, i.e., one list for
#'   each dependent variable \code{k}. Each matrix contains the hypothesized
#'   hyper-prior means of the lognormal emission distributions in each of the
#'   states in the natural (positive real numbers) scale. Hence, each
#'   matrix consists of one row (when not including covariates in the model)
#'   and \code{m} columns. If covariates are used, the number of rows in each
#'   matrix in the list is equal to 1 + n_xx_emiss (i.e., the first row
#'   corresponds to the hyper-prior means, the subsequent rows correspond to
#'   the hyper-prior values of the regression coefficients connected to each of
#'   the covariates). If covariates are used to predict the emission
#'   distribution, then \code{emiss_mu0} should be specified in the logarithmic
#'   scale, and \code{log_scale} set to \code{TRUE}.
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
#'   vector with length \code{m} containing the hypothesized variance between
#'   the subject (emission distribution) means in the natural (positive real
#'   numbers) scale, which are assumed to follow a Inverse Gamma hyper-prior
#'   distribution (note: here, the Inverse Gamma hyper-prior distribution is
#'   parametrized as a scaled inverse chi-squared distribution). If
#'   covariates are used to predict the emission distribution, then
#'   \code{emiss_V} should be specified in the logarithmic scale, and
#'   \code{log_scale} set to \code{TRUE}.
#' @param emiss_nu A list containing \code{n_dep} elements corresponding to each
#'   dependent variable \code{k}. Each element \code{k} is a numeric vector with
#'   length 1 denoting the degrees of freedom of the Inverse Gamma hyper-prior
#'   distribution on the between subject variance of the emission distribution
#'   means (note: here, the Inverse Gamma hyper-prior distribution is
#'   parametrized as a scaled inverse chi-squared distribution).
#' @param log_scale A logical scalar. Should \code{emiss_mu0} and \code{emiss_V}
#'   be specified in the logarithmic scale (\code{log_scale = TRUE}) or the
#'   natural scale (\code{log_scale = FALSE}). The default equals
#'   \code{log_scale = FALSE}.
#'
#' @return \code{prior_emiss_count} returns an object of class \code{mHMM_prior_emiss},
#'   containing informative hyper-prior values for the continuous emission
#'   distribution(s) of the multilevel hidden Markov model. The object is
#'   specifically created and formatted for use by the function \code{mHMM},
#'   and thoroughly checked for correct input dimensions.
#'   The object contains the following components:
#'   \describe{
#'   \item{\code{gen}}{A list containing the elements \code{m}, and \code{n_dep},
#'   used for checking equivalent general model properties
#'   specified under \code{prior_emiss_count} and \code{mHMM}.}
#'   \item{\code{emiss_mu0}}{A lists containing the hypothesized
#'   hyper-prior means of the lognormal distribution on
#'   the count emission probabilities.}
#'   \item{\code{emiss_K0}}{A list containing \code{n_dep} elements denoting the
#'   number of hypothetical prior subjects on which the set of hyper-prior means
#'   specified in \code{emiss_mu0} are based.}
#'   \item{\code{emiss_V}}{A list containing \code{n_dep} elements containing
#'   the variance of the Inverse Gamma hyper-prior distribution on the between
#'   subject variance of the emission distribution means.}
#'   \item{\code{emiss_nu}}{A list containing \code{n_dep} elements denoting the
#'   degrees of freedom of the Inverse Gamma hyper-prior distribution on the
#'   between subject variance of the emission distribution means.}
#'   \item{\code{n_xx_emiss}}{A numeric vector denoting the number of (level 2)
#'   covariates used to predict the emission distribution of each of the
#'   dependent variables. When no covariates are used, \code{n_xx_emiss} equals
#'   \code{NULL}.}
#'   }
#'
#' @seealso \code{\link{prior_gamma}} for manually specifying an informative
#'  hyper-prior on the transition probability matrix gamma, and
#'  \code{\link{mHMM}} for fitting a multilevel hidden Markov model.
#'
#' @examples
#' ###### Example using simulated data
#' # specifying general model properties:
#' m <- 3
#' n_dep <- 2
#'
#' # hypothesized hyper-prior values for the count emission distribution
#' manual_prior_emiss <- prior_emiss_count(
#'                         gen = list(m = m, n_dep = n_dep),
#'                         emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
#'                                          matrix(c(7, 8, 18), nrow = 1)),
#'                         emiss_K0 = list(1, 1),
#'                         emiss_V =  list(rep(16, m), rep(4, m)),
#'                         emiss_nu = list(0.1, 0.1))
#'
#' # to use the informative priors in a model, simulate multivariate count data
#' n_t     <- 100
#' n       <- 10
#'
#' # Specify group-level transition and emission means
#' gamma   <- matrix(c(0.8, 0.1, 0.1,
#'                     0.2, 0.7, 0.1,
#'                     0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
#' emiss_distr <- list(matrix(log(c( 50,
#'                               100,
#'                               150)), nrow = m, byrow = TRUE),
#'                     matrix(log(c(5,
#'                              10,
#'                              20)), nrow = m, byrow = TRUE))
#' # Simulate data
#' data_count <- sim_mHMM(n_t = n_t, n = n, data_distr = 'count', gen = list(m = m, n_dep = n_dep),
#'                   gamma = gamma, emiss_distr = emiss_distr, var_gamma = .1, var_emiss = c(.05, 0.01), log_scale = TRUE)
#'
#' # Specify starting values
#' start_gamma <- gamma
#' start_emiss <- list(matrix(c(50,
#'                              100,
#'                              150), nrow = m, byrow = TRUE),
#'                     matrix(c(5,
#'                              10,
#'                              20), nrow = m, byrow = TRUE))
#'
#'
#' # using the informative hyper-prior in a model
#' # Note that for reasons of running time, J is set at a ridiculous low value.
#' # One would typically use a number of iterations J of at least 1000,
#' # and a burn_in of 200.
#' out_3st_count_sim_infemiss <- mHMM(s_data = data_count$obs,
#'                     data_distr = "count",
#'                     gen = list(m = m, n_dep = n_dep),
#'                     start_val = c(list(start_gamma), start_emiss),
#'                     emiss_hyp_prior = manual_prior_emiss,
#'                     mcmc = list(J = 11, burn_in = 5))
#'
#' out_3st_count_sim_infemiss
#' summary(out_3st_count_sim_infemiss)
#'
#' @export
#'


prior_emiss_count <- function(gen, emiss_mu0, emiss_K0, emiss_V, emiss_nu, n_xx_emiss = NULL, log_scale = FALSE){
  if(sum(objects(gen) %in% "m") != 1 | sum(objects(gen) %in% "n_dep") != 1){
    stop("The input argument gen should contain the elements m and n_dep")
  } else if (is.null(log_scale)){
    stop("The input argument log_scale should be either FALSE for specification in the natural scale, or TRUE for specification in the logarithmic scale.")
  } else if(!log_scale %in% c(FALSE, TRUE)){
    stop("The input argument log_scale should be either FALSE for specification in the natural scale, or TRUE for specification in the logarithmic scale.")
  } else if(log_scale == FALSE & !is.null(n_xx_emiss)){
    stop("According to n_xx_emiss covariates have been used to predict the count emission distribution. Please use the logarithmic scale to specify emiss_mu0 and emiss_V, and set `log_scale = TRUE`.")
  }
  m <- gen$m
  n_dep <- gen$n_dep
  if(is.null(n_xx_emiss)){
    n_xx_int <- rep(1, n_dep)
  } else {
    if(length(n_xx_emiss) != n_dep){
      stop(paste("n_xx_emiss should be a numeric vector with length n_dep, here", n_dep, "."))
    }
    n_xx_int <- n_xx_emiss + 1
  }

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
  for(k in 1:n_dep){
    if(sum(m == sapply(emiss_mu0, dim)[2,]) != n_dep){
      stop(paste("The matrix relating to dependent variable", k, "of the input argument emiss_mu0 should consist of m, here", m, ", columns."))
    }
    if(n_xx_int[k] == 1 & dim(emiss_mu0[[k]])[1] != 1){
      stop(paste("According to the input paramter n_xx_emiss no covariates are used to predict the emission distribution of dependent variable", k, ". Hence, within input argument emiss_mu0, the matrix relating to dependent variable", k, ", should contain 1 row."))
    }

    if(n_xx_int[k] > 1 & dim(emiss_mu0[[k]])[1] != n_xx_int[k]){
      stop(paste("According to the input paramter n_xx_emiss", n_xx_emiss[k], "covariates are used to predict the emission distribution of dependent variable", k, ". Hence, within input argument emiss_mu0, the matrix relating to dependent variable", k, ", should contain 1 + n_xx_emiss =", 1 + n_xx_emiss[k], "rows."))
    }
  }

  #### checking emiss_K0 ####
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

  #### put parameters in the correct scale ####
  if(log_scale == FALSE){
    emiss_V <- var_to_logvar(gen = gen, emiss_mu = emiss_mu0, var_emiss = emiss_V, byrow = TRUE)
    for(q in 1:n_dep){
      # emiss_V[[q]]     <- var_to_logvar(emiss_mu0[[q]][1,], emiss_V[[q]])
      emiss_mu0[[q]]   <- log(emiss_mu0[[q]])
    }
  }

  #### return output ####
  out <- list(gen = gen, emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0, emiss_nu = emiss_nu, emiss_V = emiss_V, n_xx_emiss = n_xx_emiss)
  class(out) <- append(class(out), c("mHMM_prior_emiss", "count"))
  return(out)
}
