#' Transform the between-subject variance in the positive scale to the logvarince in the logarithmic scale
#'
#' \code{var_to_logvar} returns the desired between-subject logvariance in the
#'   logarithmic scale corresponding to the between-subject variance in the
#'   (real) positive scale specified as input. The logvariance is used as the
#'   dispersion parameter in the lognormal distribution adopted as prior for
#'   the subject-specific Poisson means in the implementation of \code{mHMM}
#'   and \code{\link{sim_mHMM}} of \code{\link{mHMMbayes}} for count
#'   distributed data. It takes as inputs the group-level Poisson means
#'   \code{emiss_mu} and the desired between-subject variance \code{var_emiss},
#'   both in the (real) positive scale.
#'
#' @param gen List containing the following elements denoting the general model
#'   properties:
#'   \itemize{
#'   \item{\code{m}: numeric vector with length 1 denoting the number
#'   of hidden states}
#'   \item{\code{n_dep}: numeric vector with length 1 denoting the
#'   number of dependent variables}}
#' @param emiss_mu A list containing \code{n_dep} matrices, i.e., one list for
#'   each dependent variable \code{k}. Each matrix contains the group-level
#'   means of the lognormal emission distributions in each of the
#'   states in the natural (positive real numbers) scale. Hence, each
#'   matrix consists of either 1 row (when not including covariates in the
#'   model) and \code{m} columns (see \code{emiss_mu0} in
#'   \code{\link{prior_emiss_count}}), or \code{m} rows and 1 column
#'   (when not including covariates in the model; see \code{emiss_distr} in
#'   \code{\link{sim_mHMM}}), denoting the mean of state
#'   \emph{i} (column or row \emph{i}) of the lognormal distribution used as
#'   prior for the Poisson emissions. By default it is assumed that the
#'   matrices contain 1 row and \code{m} columns, as specified for
#'   \code{emiss_mu0} in \code{\link{prior_emiss_count}} (see argument
#'   \code{byrow}). If \code{emiss_mu} were to be specified using \code{m} rows
#'   and 1 column as for \code{emiss_distr} in \code{\link{sim_mHMM}}, then set
#'   the argument \code{byrow = FALSE}.
#'
#'   Note: in every case the means should be specified in the natural
#'   (real positive) scale.
#' @param var_emiss A list containing \code{n_dep} elements corresponding to
#'   each of the dependent variables \code{k}, where each element \code{k} is a
#'   vector with length \code{m} denoting the amount of variance between the
#'   subject (emission distribution) means in the natural (positive real
#'   numbers) scale. It follows a similar specification as \code{emiss_V} in
#'   \code{\link{prior_emiss_count}}.
#' @param byrow A logical scalar indicating whether the emission means are
#'   specified in the first row (\code{byrow = TRUE}) or the first column
#'   (\code{byrow = FALSE}) of the wah of the \code{n_dep} matrices listed in
#'   \code{emiss_mu}. Use \code{byrow = TRUE} if \code{emiss_mu} is entered
#'   as \code{emiss_mu0} in \code{\link{prior_emiss_count}}, and
#'   \code{byrow = FALSE} if \code{emiss_mu} is entered as \code{emiss_distr}
#'   in \code{\link{sim_mHMM}}. By default, \code{byrow = TRUE}.
#'
#' @return \code{var_to_logvar} Returns a list of \code{n_dep} numeric vectors
#'   of \code{m} elements denoting the state \code{i}-specific logvariance
#'   (between-subject variance in the logarithmic scale) for the
#'   \code{k}-dependent variable used as dispersion parameter in the lognormal
#'   prior for the Poisson emission distribution.
#'
#' @seealso \code{\link{mHMM}} for fitting the multilevel hidden Markov model,
#'   and \code{\link{sim_mHMM}} for simulating data from a multilevel hidden
#'   Markov model.
#'
#'
#'
#' @examples
#' ###### Example on package data
#' \donttest{
#'
#' ## Example: specifying count priors manually using both positive and log scale:
#'
#' # Define general parameters:
#' m       <- 3        # Number of hidden states
#' n_dep   <- 2        # Number of dependent variables
#'
#' # Specify priors manually on the positive scale for emuss_mu0 and emiss_V:
#' manual_prior_emiss_1 <- prior_emiss_count(
#'   gen = list(m = m, n_dep = n_dep),
#'   emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
#'                    matrix(c(7, 8, 18), nrow = 1)),
#'   emiss_K0 = list(1, 1),
#'   emiss_V =  list(c(16,25,32),
#'                   rep(4, m)),
#'   emiss_nu = list(0.1, 0.1),
#'   log_scale = FALSE)
#'
#' # Define logmu and logvar:
#' logmu <-  list(matrix(log(c(30, 70, 170)), nrow = 1),
#'                matrix(log(c(7, 8, 18)), nrow = 1))
#'
#' logvar <- var_to_logvar(gen = list(m = m, n_dep = n_dep),
#'                         emiss_mu = list(matrix(c(30, 70, 170), nrow = 1),
#'                                         matrix(c(7, 8, 18), nrow = 1)),
#'                         var_emiss = list(c(16,25,32),
#'                                          rep(4, m)),
#'                         byrow = TRUE)
#'
#' # Specify priors manually on the log scale for emiss_mu0 and emiss_V:
#' manual_prior_emiss_2 <- prior_emiss_count(
#'   gen = list(m = m, n_dep = n_dep),
#'   emiss_mu0 = logmu,
#'   emiss_K0 = list(1, 1),
#'   emiss_V =  logvar,
#'   emiss_nu = list(0.1, 0.1),
#'   log_scale = TRUE)
#'
#' # Check whether they are identical:
#' identical(manual_prior_emiss_1, manual_prior_emiss_2)
#'
#' }
#'
#' @export

var_to_logvar <- function(gen, emiss_mu, var_emiss, byrow = TRUE){

  # Initialize parameters
  n_dep   <- gen$n_dep
  m       <- gen$m
  logvar  <- vector("list", n_dep)

  # Insert checks (stops/warnings)
  if(sum(objects(gen) %in% "m") != 1 | sum(objects(gen) %in% "n_dep") != 1){
    stop("The input argument gen should contain the elements m and n_dep")
  }

  #### checking emiss_mu ####
  if (any(unlist(emiss_mu) <= 0)){
    stop(paste("All values in emiss_mu should be real positive numbers (emiss_mu > 0)."))
  }
  if(any(!sapply(emiss_mu, is.matrix))){
    stop(paste("emiss_mu should be a list containing", n_dep, "matrices; one matrix for each dependent variable."))
  }
  if(byrow){

    # State means in 1st row
    if(!is.list(emiss_mu)){
      stop(paste("emiss_mu should be a list containing", n_dep, "matrices; one matrix for each dependent variable."))
    }
    if(length(emiss_mu) != n_dep ){
      stop(paste("emiss_mu should be a list containing", n_dep, "matrices; one matrix for each dependent variable."))
    }
    for(q in 1:n_dep){
      if (dim(emiss_mu[[q]])[2] != m){
        stop(paste("The number of columns of the emission distribution matrix in element", q, "should be
             equal to the number of states, which is", m, "."))
      }
    }

  } else {

    # State means in 1st column
    if(!is.list(emiss_mu)){
      stop("The format of emiss_mu should be a list with", n_dep, "elements.")
    }
    if(length(emiss_mu) != n_dep){
      stop("The number of dependent variables specified in n_dep and the number of elements specified in the list emiss_mu should be equal")
    }
    for(q in 1:n_dep){
      if (dim(emiss_mu[[q]])[1] != m){
        stop(paste("The number of rows of the emission distribution matrix in element", q, "should be
             equal to the number of states, which is", m, "."))
      }
    }

  }

  #### checking emiss_V ####
  if(!is.list(var_emiss)){
    stop(paste("var_emiss should be a list containing n_dep, here", n_dep,", elements, where each element is vector with lenght m, here", m, "."))
  }
  if(length(var_emiss) != n_dep | sum(m == sapply(var_emiss, length)) != n_dep){
    stop(paste("var_emiss should be a list containing n_dep, here", n_dep,", elements, where each element is vector with lenght m, here", m, "."))
  }
  if (any(unlist(var_emiss) <= 0)){
    stop(paste("All values in var_emiss should be real non-negative numbers (var_emiss larger or equal to 0."))
  }

  # Cycle over dependent variables
  for(q in 1:n_dep){
    if(byrow){
      logmu <- log(as.numeric(emiss_mu[[q]][1,]))
    } else {
      logmu <- log(as.numeric(emiss_mu[[q]][,1]))
    }
    logvar[[q]] <- log(0.5*exp(-2*logmu)*(exp(2*logmu) + sqrt(4*exp(2*logmu)*var_emiss[[q]]+exp(4*logmu))))
  }

  return(logvar)

}
