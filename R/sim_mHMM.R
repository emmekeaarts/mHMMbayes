#' Simulate data using a multilevel hidden Markov model
#'
#' \code{sim_mHMM} simulates data for multiple subjects, for which the data have
#' either categorical or continuous (i.e., normally distributed)
#' observations that follow a hidden Markov model (HMM) with a
#' multilevel structure. The multilevel structure implies that each subject is
#' allowed to have its own set of parameters, and that the parameters at the
#' subject level (level 1) are tied together by a population distribution at
#' level 2 for each of the corresponding parameters. The shape of the population
#' distribution for each of the parameters is a normal distribution. In addition
#' to (natural and/or unexplained) heterogeneity between subjects, the subjects
#' parameters can also depend on a covariate.
#'
#' In simulating the data, having a multilevel structure means that the
#' parameters for each subject are sampled from the population level
#' distribution of the corresponding parameter. The user specifies the
#' population distribution for each parameter: the average population transition
#' probability matrix and its variance, and the average population emission
#' distribution and its variance. For now, the variance of the mean population
#' parameters is assumed fixed for all components of the transition probability
#' matrix and for all components of the emission distribution.
#'
#' One can simulate multivariate data. That is, the hidden states depend on more
#' than 1 observed variable simultaneously. The distributions of multiple
#' dependent variables for multivariate data are assumed to be independent, and
#' all distributions for one dataset have to be of the same type (either
#' categorical or continuous).
#'
#' Note: the subject specific) initial state distributions (i.e., the
#' probability of each of the states at the first time point) needed to simulate
#' the data are obtained from the stationary distributions of the subject
#' specific transition probability matrices gamma.
#'
#'
#' \code{beta}: As the first element in each row of \code{gamma} is used as
#' reference category in the Multinomial logistic regression, the first matrix
#' in the list \code{beta} used to predict transition probability matrix
#' \code{gamma} has a number of rows equal to \code{m} and the number of columns
#' equal to \code{m} - 1. The first element in the first row corresponds to the
#' probability of switching from state one to state two. The second element in
#' the first row corresponds to the probability of switching from state one to
#' state three, and so on. The last element in the first row corresponds to the
#' probability of switching from state one to the last state. The same principle
#' holds for the second matrix in the list \code{beta} used to predict
#' categorical emission distribution(s) \code{emiss_distr}: the first element in
#' the first row corresponds to the probability of observing category two in
#' state one. The second element in the first row corresponds to the probability
#' of observing category three is state one, and so on. The last element in the
#' first row corresponds to the probability of observing the last category in
#' state one.
#'
#' @inheritParams mHMM
#' @param n_t Numeric vector with length 1 denoting the length of the observed
#'   sequence to be simulated for each subject. To only simulate subject
#'   specific transition probability matrices gamma and emission distributions
#'   (and no data), set \code{t} to 0.
#' @param n Numeric vector with length 1 denoting the number of subjects for
#'   which data is simulated.
#' @param data_distr String vector with length 1 denoting the observation type
#'   of the data to be simulated. Should be set to either \code{'categorical'}
#'   or \code{'continuous'}. Note that when simulating multivariate data, all
#'   dependent variables are assumed to be of the same observation type. The
#'   default equals to \code{data_distr = 'categorical'}.
#' @param m The argument \code{m} is deprecated; please specify using the input
#'   parameter \code{gen}.
#' @param n_dep The argument \code{n_dep} is deprecated; please specify using
#'   the input parameter \code{gen}.
#' @param q_emiss The argument \code{q_emiss} is deprecated; please specify
#'   using the input parameter \code{gen} (only to be specified when simulating
#'   categorical data).
#' @param start_state Optional numeric vector with length 1 denoting in which
#'   state the simulated state sequence should start. If left unspecified, the
#'   simulated state for time point 1 is sampled from the initial state
#'   distribution (which is derived from the transition probability matrix
#'   gamma).

#' @param gamma A matrix with \code{m} rows and \code{m} columns containing the
#'   average population transition probability matrix used for simulating the
#'   data. That is, the probability to switch from hidden state \emph{i} (row
#'   \emph{i}) to hidden state \emph{j} (column  \emph{j}).
#' @param emiss_distr A list with \code{n_dep} elements containing the average
#'   population emission distribution(s) of the observations given the hidden
#'   states for each of the dependent variables. If \code{data_distr =
#'   'categorical'}, each element is a matrix with \code{m} rows and
#'   \code{q_emiss[k]} columns for each of the \code{k} in \code{n_dep} emission
#'   distribution(s). That is, the probability of observing category \emph{q}
#'   (column \emph{q}) in state \emph{i} (row \emph{i}). If \code{data_distr =
#'   'continuous'}, each element is a matrix with \code{m} rows and 2 columns;
#'   the first column denoting the mean of state \emph{i} (row \emph{i}) and the
#'   second column denoting the standard deviation of state \emph{i}
#'   (row \emph{i}) of the Normal distribution.
#' @param xx_vec List of 1 + \code{n_dep} vectors containing the covariate(s) to
#'   predict the transition probability matrix \code{gamma} and/or (specific)
#'   emission distribution(s) \code{emiss_distr} using the regression parameters
#'   specified in \code{beta} (see below). The first element in the list
#'   \code{xx_vec} is used to predict the transition matrix. Subsequent elements
#'   in the list are used to predict the emission distribution of (each of) the
#'   dependent variable(s). This means that the covariate used to predict
#'   \code{gamma} and \code{emiss_distr} can either be the same covariate,
#'   different covariates, or a covariate for certain elements and none for the
#'   other. At this point, it is only possible to use one covariate for both
#'   \code{gamma} and \code{emiss_distr}. For all elements in the list,
#'   the number of observations in the vectors should be  equal to the number of
#'   subjects to be simulated \code{n}. If \code{xx_vec} is omitted completely,
#'   \code{xx_vec} defaults to NULL, resembling no covariates at all. Specific
#'   elements in the list can also be left empty (i.e., set to \code{NULL}) to
#'   signify that either the transition probability matrix or (one of) the
#'   emission distribution(s) is not predicted by covariates.
#' @param beta List of 1 + \code{n_dep} matrices containing the regression
#'   parameters to predict \code{gamma} and/or \code{emiss_distr} in combination
#'   with \code{xx_vec} using (Multinomial logistic) regression. The first
#'   matrix is used to predict the transition probability matrix \code{gamma}.
#'   The subsequent matrices are used to predict the emission distribution(s)
#'   \code{emiss_distr} of the dependent variable(s). For \code{gamma} and
#'   categorical emission distributions, one regression parameter is specified
#'   for each element in \code{gamma} and \code{emiss_distr}, with the following
#'   exception. The first element in each row of \code{gamma} and/or
#'   \code{emiss_distr} is used as reference category in the Multinomial
#'   logistic regression. As such, no regression parameters can be specified for
#'   these parameters. Hence, the first element in the list \code{beta} to
#'   predict \code{gamma} consist of a matrix with the number of rows equal to
#'   \code{m} and the number of columns equal to \code{m} - 1. For categorical
#'   emission distributions, the subsequent elements in the list \code{beta} to
#'   predict \code{emiss_distr} consist of a matrix with the number of rows
#'   equal to \code{m} and the number of columns equal to \code{q_emiss[k]} - 1
#'   for each of the \code{k} in \code{n_dep} emission distribution(s). See
#'   \emph{details} for more information. For continuous emission distributions,
#'   the subsequent elements in the list \code{beta} consist of a matrix with
#'   the number of rows equal to \code{m} and 1 column.
#'
#'   Note that if \code{beta} is specified, \code{xx_vec} has to be specified as
#'   well. If \code{beta} is omitted completely, \code{beta} defaults to NULL,
#'   resembling no prediction of \code{gamma} and \code{emiss_distr} using
#'   covariates. One of the elements in the list can also be left empty
#'   (i.e., set to \code{NULL}) to signify that either the transition
#'   probability matrix or a specific emission distribution is not predicted by
#'   covariates.
#' @param var_gamma A numeric vector with length 1 denoting the amount of variance between subjects in
#'   the transition probability matrix. Note that this value corresponds to the
#'   variance of the parameters of the Multinomial distribution (i.e., the
#'   intercepts of the regression equation of the Multinomial distribution used
#'   to sample the transition probability matrix), see details below. In
#'   addition, only one variance value can be specified for the complete
#'   transition probability matrix, hence the variance is assumed fixed across
#'   all components. The default equals 0.1, which corresponds to little
#'   variation between subjects. If one wants to simulate data from exactly the
#'   same HMM for all subjects, var_gamma should be set to 0. Note that if data
#'   for only 1 subject is simulated (i.e., n = 1), \code{var_gamma} is set to
#'   0.
#' @param var_emiss A numeric vector with length \code{n_dep} denoting the
#'   amount of variance between subjects in the emission distribution(s). For
#'   categorical data, this value corresponds to the variance of the parameters
#'   of the Multinomial distribution (i.e., the intercepts of the regression
#'   equation of the Multinomial distribution used to sample the components of
#'   the emission distribution), see details below. For continuous data, this
#'   value corresponds to the variance in the mean of the emission
#'   distribution(s) across subjects. Note that only one variance value can be
#'   specified each emission distribution, hence the variance is assumed fixed
#'   across states (and, for the categorical distribution, categories within
#'   a state) within an emission distribuitonß. The default equals 0.1,
#'   which corresponds to little variation between subjects given categorical
#'   observations. If one wants to simulate data from exactly the same HMM for
#'   all subjects, var_emiss should be set to a vector of 0's. Note that if data
#'   for only 1 subject is simulated (i.e., n = 1), \code{var_emiss} is set to a
#'   vector of 0's.
#' @param return_ind_par A logical scalar. Should the subject specific
#'   transition probability matrix \code{gamma} and emission probability matrix
#'   \code{emiss_distr} be returned by the function (\code{return_ind_par =
#'   TRUE}) or not (\code{return_ind_par = FALSE}). The default equals
#'   \code{return_ind_par = FALSE}.
#'
#'
#' @return The following components are returned by the function \code{sim_mHMM}:
#' \describe{
#'   \item{\code{states}}{A matrix containing the simulated hidden state
#'   sequences, with one row per hidden state per subject. The first column
#'   indicates subject id number. The second column contains the simulated
#'   hidden state sequence, consecutively for all subjects. Hence, the id number
#'   is repeated over the rows (with the number of repeats equal to the length
#'   of the simulated hidden state sequence \code{T} for each subject).}
#'   \item{\code{obs}}{A matrix containing the simulated observed outputs, with
#'   one row per simulated observation per subject. The first column indicates
#'   subject id number. The second column contains the simulated observation
#'   sequence, consecutively for all subjects. Hence, the id number is repeated
#'   over rows (with the number of repeats equal to the length of the simulated
#'   observation sequence \code{T} for each subject).}
#'   \item{\code{gamma}}{A list containing \code{n} elements with the simulated
#'   subject specific transition probability matrices \code{gamma}. Only
#'   returned if \code{return_ind_par} is set to \code{TRUE}.}
#'   \item{\code{emiss_distr}}{A list containing \code{n} elements with the
#'   simulated subject specific emission probability matrices
#'   \code{emiss_distr}. Only returned if \code{return_ind_par} is set to
#'   \code{TRUE}.}
#' }
#'
#'
#' @seealso \code{\link{mHMM}} for analyzing multilevel hidden Markov data.
#'
#'
#' @examples
#' # simulating data for 10 subjects with each 100 observations
#' n_t     <- 100
#' n       <- 10
#' m       <- 3
#' n_dep   <- 1
#' q_emiss <- 4
#' gamma   <- matrix(c(0.8, 0.1, 0.1,
#'                     0.2, 0.7, 0.1,
#'                     0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
#' emiss_distr <- list(matrix(c(0.5, 0.5, 0.0, 0.0,
#'                              0.1, 0.1, 0.8, 0.0,
#'                              0.0, 0.0, 0.1, 0.9), nrow = m, ncol = q_emiss, byrow = TRUE))
#' data1 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                   gamma = gamma, emiss_distr = emiss_distr, var_gamma = 1, var_emiss = 1)
#' head(data1$obs)
#' head(data1$states)
#'
#' # including a covariate to predict (only) the transition probability matrix gamma
#' beta      <- rep(list(NULL), 2)
#' beta[[1]] <- matrix(c(0.5, 1.0,
#'                      -0.5, 0.5,
#'                       0.0, 1.0), byrow = TRUE, ncol = 2)
#' xx_vec      <- rep(list(NULL),2)
#' xx_vec[[1]] <-  c(rep(0,5), rep(1,5))
#' data2 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                   gamma = gamma, emiss_distr = emiss_distr, beta = beta, xx_vec = xx_vec,
#'                   var_gamma = 1, var_emiss = 1)
#'
#'
#' # simulating subject specific transition probability matrices and emission distributions only
#' n_t <- 0
#' n <- 5
#' m <- 3
#' n_dep   <- 1
#' q_emiss <- 4
#' gamma <- matrix(c(0.8, 0.1, 0.1,
#'                   0.2, 0.7, 0.1,
#'                   0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
#' emiss_distr <- list(matrix(c(0.5, 0.5, 0.0, 0.0,
#'                              0.1, 0.1, 0.8, 0.0,
#'                              0.0, 0.0, 0.1, 0.9), nrow = m, ncol = q_emiss, byrow = TRUE))
#' data3 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                   gamma = gamma, emiss_distr = emiss_distr, var_gamma = 1, var_emiss = 1)
#' data3
#'
#' data4 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                   gamma = gamma, emiss_distr = emiss_distr, var_gamma = .5, var_emiss = .5)
#' data4
#'
#' # simulating multivariate continuous data
#' n_t     <- 100
#' n       <- 10
#' m       <- 3
#' n_dep   <- 2
#'
#' gamma   <- matrix(c(0.8, 0.1, 0.1,
#'                     0.2, 0.7, 0.1,
#'                     0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
#'
#' emiss_distr <- list(matrix(c( 5, 1,
#'                              10, 1,
#'                              15, 1), nrow = m, byrow = TRUE),
#'                      matrix(c(0.5, 0.2,
#'                               1.0, 0.5,
#'                               2.0, 0.3), nrow = m, byrow = TRUE))
#'
#' data_cont <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep),
#'                       data_distr = 'continuous', gamma = gamma, emiss_distr = emiss_distr,
#'                       var_gamma = .5, var_emiss = c(.5, 0.01))
#'
#' head(data_cont$states)
#' head(data_cont$obs)


#' @export

sim_mHMM <- function(n_t, n, data_distr = 'categorical', gen, gamma, emiss_distr, start_state = NULL,
                     xx_vec = NULL, beta = NULL, var_gamma = 0.1, var_emiss = NULL,
                     return_ind_par = FALSE, m, n_dep, q_emiss){

  if(!missing(m)){
    warning("The argument m is deprecated; please specify using the input parameter gen.")
  }
  if(!missing(n_dep)){
    warning("The argument n_dep is deprecated; please specify using the input parameter gen.")
  }
  if(!missing(q_emiss)){
    warning("The argument q_emiss is deprecated; please specify using the input parameter gen.")
  }

  if(!missing(gen)){
    if(sum(objects(gen) %in% "m") != 1 | sum(objects(gen) %in% "n_dep") != 1 | (sum(objects(gen) %in% "q_emiss") != 1 & data_distr == 'categorical')){
      stop("The input argument gen should contain the elements m, n_dep and q_emiss.")
    }
    m <- gen$m
    n_dep <- gen$n_dep
    if(data_distr == 'categorical'){
      q_emiss <- gen$q_emiss
    }
  }

  if(missing(m) & missing(gen)){
    stop("Please specify the number of hidden states m via the input parameter gen.")
  }
  if(missing(n_dep) & missing(gen)){
    n_dep <- 1
    warning("Please specify the number of dependent variables n_dep via the input parameter gen. Model now assums n_dep = 1")
  }
  if(data_distr == 'categorical'){
    if(missing(q_emiss) & missing(gen)){
      stop("Please specify the number of observed categories for each categorical emission distribtution q_emiss via the input parameter gen.")
    }
  }

  if (dim(gamma)[1] != m | dim(gamma)[2] != m){
    stop(paste("The transiton probability matrix gamma should be a", m, "by", m, "matrix."))
  }
  if(!isTRUE(all.equal(apply(gamma,1,sum), rep(1,m)))){
    stop("The elements in each row of the transition probability matrix gamma should sum up to 1")
  }
  if(!is.list(emiss_distr)){
    stop("The format of emiss_distr should be a list with", n_dep, "elements.")
  }
  if(length(emiss_distr) != n_dep){
    stop("The number of dependent variables specified in n_dep and the number of elements specified in the list emiss_distr should be equal")
  }
  if(data_distr == 'categorical'){
    if(length(q_emiss) != n_dep){
      stop("The lenght of q_emiss specifying the number of output categories for each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
    }
  }
  for(i in 1:n_dep){
    if (dim(emiss_distr[[i]])[1] != m){
      stop(paste("The number of rows of emission distribution matrix in element", i, "should be
             equal to the number of states, which is", m, "."))
    }
    if(data_distr == 'categorical'){
     if (dim(emiss_distr[[i]])[2] != q_emiss[i]){
       stop(paste("The number of columns of the emission distribution matrix should be
                 equal to the number of observable categories, which is", q_emiss[i], ". See emission distribution in element", i, "."))
      }
      if(!isTRUE(all.equal(apply(emiss_distr[[i]], 1, sum), rep(1, m)))){
       stop("The elements in each row of the emission distribution matrix should sum up to 1, see emission distribution in element", i, ".")
      }
    }
#   if(data_distr == 'continuous'){
#
#   }
  }
  if((is.null(xx_vec) & !is.null(beta)) | (!is.null(xx_vec) & is.null(beta))){
    stop("Either only xx_vec or only beta is specified. Please specify both 1) the values for the covariate
         in xx_vec and 2) the values of the regression parameters in beta, to allow correct simulation of the
         data.")
  }
  if(!is.null(xx_vec)){
    if((!is.null(xx_vec[[1]]) & is.null(beta[[1]])) |
       (!is.null(xx_vec[[2]]) & is.null(beta[[2]]))){
      stop("Either only xx_vec or only beta is specified in one of the elements.
            Please specify both 1) the values for the covariate in xx_vec and 2)
            the values of the regression parameters in beta if either one is not
            empty, to allow correct simulation of the data.")
    }
  }
  if(!is.null(beta)){
  # extend to all 1 + n_dep
      if((!is.null(beta[[1]]) & is.null(xx_vec[[1]])) |
       (!is.null(beta[[2]]) & is.null(xx_vec[[2]]))){
      stop("Either only xx_vec or only beta is specified in one of the elements.
            Please specify both 1) the values for the covariate in xx_vec and 2)
            the values of the regression parameters in beta if either one is not
            empty, to allow correct simulation of the data.")
    }
  }
  if(!is.null(xx_vec)){
  # extend to all 1 + n_dep
    if((!is.null(xx_vec[[1]]) & length(xx_vec[[1]]) != n) |
      (!is.null(xx_vec[[2]]) & length(xx_vec[[2]]) != n)){
      stop("The length of the vectors in xx_vec should be equal to the number of subjects to be simulated,
         set in n, if (the element in) xx_vec is not set to NULL.")
    }
  }
  if (!is.null(beta)){
    if (!is.null(beta[[1]])){
      if ((dim(beta[[1]])[1] != (m)) | (dim(beta[[1]])[2] != (m-1))){
      stop(paste("The first element of beta to predict the transiton probability matrix gamma should be a m (", m, " ) by m - 1 (", m - 1, ") matrix."))
      }
    }
    if (!is.null(beta[[2]]) & data_distr == 'categorical'){
    # extend to all 1 + n_dep and continuous
      if((dim(beta[[2]])[1] != (m)) | (dim(beta[[2]])[2] != (q_emiss[1]-1))){
      stop(paste("The second element of beta to predict the emission distribution should be a m (", m, ") by q_emiss - 1 (", q_emiss[1] - 1, ") matrix."))
      }
    }
  }
  if(is.null(xx_vec)){
    xx_vec <- rep(list(NULL), n_dep + 1)
    for(i in 1:(n_dep + 1)){
      xx_vec[[i]] <- rep(1,n)
    }
  } else {
    for(i in 1:(n_dep + 1)){
      if(is.null(xx_vec[[i]])) {
        xx_vec[[i]] <- rep(1,n)
      }
    }
  }
  if(is.null(beta)){
    beta <- rep(list(NULL), n_dep + 1)
    beta[[1]] <- matrix(0, ncol = m - 1, nrow = m)
    for(i in 2:(n_dep + 1)){
      if(data_distr == 'categorical'){
        beta[[i]] <- matrix(0, ncol = q_emiss[i-1] - 1, nrow = m)
      } else if (data_distr == 'continuous'){
        beta[[i]] <- matrix(0, ncol = 1, nrow = m)
      }
    }
  } else {
    if(is.null(beta[[1]])) {
      beta[[1]] <- matrix(0, ncol = m - 1, nrow = m)
    }
    for (i in 2:(n_dep + 1)){
      if (is.null(beta[[i]])) {
        if(data_distr == 'categorical'){
          beta[[i]] <- matrix(0, ncol = q_emiss[i-1] - 1, nrow = m)
        } else if (data_distr == 'continuous'){
          beta[[i]] <- matrix(0, ncol = 1, nrow = m)
        }
      }
    }
  }

  if(n == 1){
    var_gamma <- 0
    var_emiss <- rep(0, n_dep)
  }
  if(is.null(var_emiss)){
    var_emiss <- rep(0.1, n_dep)
  } else if(length(var_emiss) != n_dep){
    stop("The lenght of var_emiss specifying variance between subjects in each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
  }

  #############
  # Simulating the data ---------------------
  #############

  states <- matrix(ncol = 2, nrow = n_t*n)
  states[,1] <- rep(1:n, each = n_t)
  obs <- matrix(ncol = 1 + n_dep, nrow = n_t*n)
  obs[,1] <- rep(1:n, each = n_t)
  sub_gamma <- rep(list(NULL), n)
  sub_emiss <- rep(list(vector("list", n_dep)), n)
  mnl_gamma <- prob_to_int(gamma)
  if(data_distr == "categorical"){
    mnl_emiss <- rep(list(NULL), n_dep)
    for(i in 1:n_dep){
      mnl_emiss[[i]] <- prob_to_int(emiss_distr[[i]])
    }
  }
  for(j in 1:n){
    sub_gamma[[j]] <- int_to_prob(mnl_gamma + xx_vec[[1]][j] * beta[[1]] +
                                    rnorm(n = m * (m-1), mean = 0, sd = sqrt(var_gamma)))
    for(i in 1:n_dep){
      if(data_distr == "categorical"){
        sub_emiss[[j]][[i]] <- int_to_prob(mnl_emiss[[i]] + xx_vec[[1+i]][j] * beta[[1+i]] +
                                           rnorm(n = m * (q_emiss[i]-1), mean = 0, sd = sqrt(var_emiss[i])))
      } else if(data_distr == "continuous"){
        sub_emiss[[j]][[i]] <- emiss_distr[[i]]
        sub_emiss[[j]][[i]][,1] <- emiss_distr[[i]][,1] +  xx_vec[[1+i]][j] * beta[[1+i]] +
        rnorm(n = m, mean = 0, sd = sqrt(var_emiss[i]))
      }
    }

    if(n_t != 0){
      init <- solve(t(diag(m) - sub_gamma[[j]] + 1), rep(1, m))
      if (is.null(start_state)){
        states[((j-1) * n_t + 1), 2] <- sample(x = 1:m, size = 1, prob = init)
      } else {
        states[((j-1) * n_t + 1), 2] <- start_state
      }
      if(data_distr == "categorical"){
        for(i in 1:n_dep){
          obs[((j-1) * n_t + 1), (1+i)] <- sample(x = 1:q_emiss[i], size = 1, prob = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],])
        }
      } else if (data_distr == "continuous"){
        for(i in 1:n_dep){
          obs[((j-1) * n_t + 1), (1+i)] <- rnorm(1, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],1],
                                                 sd = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],2])
        }
      }
      for(t in 2:n_t){
        states[((j-1) * n_t + t), 2] <- sample(x = 1:m, size = 1, prob = sub_gamma[[j]][states[((j-1) * n_t + t - 1), 2],])
        if(data_distr == "categorical"){
          for(i in 1:n_dep){
            obs[((j-1) * n_t + t), (1+i)] <- sample(x = 1:q_emiss[i], size = 1, prob = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],])
          }
        } else if (data_distr == "continuous"){
          for(i in 1:n_dep){
            obs[((j-1) * n_t + t), (1+i)] <- rnorm(1, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],1],
                                                   sd = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],2])
          }
        }
      }
    }
  }

  #############
  # Returning output  ---------------------
  #############
  colnames(states) <- c("subj", "state")
  colnames(obs)    <- c("subj", paste("observation", 1:n_dep))
  if (return_ind_par == FALSE & n_t != 0){
    return(list(states = states, obs = obs))
  } else if (return_ind_par == TRUE & n_t != 0){
    return(list(states = states, obs = obs, subject_gamma = sub_gamma, subject_emiss = sub_emiss))
  } else if (n_t == 0){
    return(list(subject_gamma = sub_gamma, subject_emiss = sub_emiss))
  }
}
