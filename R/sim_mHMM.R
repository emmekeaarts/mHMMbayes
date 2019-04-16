#' Simulate data using a multilevel hidden markov model
#'
#' This function simulates data for multiple subjects, for which the data have
#' categorical observations that follow a hidden Markov model (HMM) with an
#' multilevel structure. The multilevel structure implies that each subject is
#' allowed to have it's own set of parameters, and that the parameters at the
#' subject level (level 1) are tied together by a population distribution at
#' level 2 for each of the corresponding parameters. The shape of the population
#' distribution for each of the parameters is a normal (i.e., Gaussian)
#' distribution. In addition to (natural and/or unexplained) heterogeneity
#' between subjects, the subjects parameters can also depend on a (set of)
#' covariate(s).
#'
#' In simulating the data, having a multilevel structure means that the
#' parameters for each subject are sampled from the population level
#' distribution of the corresponding parameter. The user specifies the
#' population distribution for each parameter: the average population transition
#' probability matrix and its variance, and the average population emission
#' distribution and its variance. For now, the variance is assumed fixed for all
#' components of the transition probability matrix and for all components of the
#' emissiondistribution, and the simulated data can only consist of one
#' dependent variable. In additon, at this point only one dependent variable can
#' be simulated. That is, the hidden Markov model is a univariate hidden Markov
#' model.
#'
#' Note: the subject specific) initial state distributions (i.e., the
#' probability of each of the states at the first time point) needed to simulate
#' the data are obtained from the stationary distributions of the subject
#' specific transition probability matrices gamma.
#'
#' \subsection{Covariates}{ Here is some info on using covariates ... (to be
#' completed)
#'
#' \code{beta}: The first element in the list \code{beta} to predict
#' \code{gamma} consist of a matrix with the number of rows equal to \code{m}
#' and the number of columns equal to \code{m} - 1. The first element in the
#' first row corresponds to the probability of switching from state one to state
#' two. The second element in the first row corresponds to the probability of
#' switching from state one to state three, and so on. The last element in the
#' first row corresponds to the probability of switching from state one to the
#' last state. The second element in the list \code{beta} to predict
#' \code{emiss_distr} consist of a matrix with the number of rows equal to
#' \code{m} and the number of columns equal to \code{pr} - 1. The first element
#' in the first row corresponds to the probability of observing category two in
#' state one. The second element in the first row corresponds to the probability
#' of observing category three is state one, and so on. The last element in the
#' first row corresponds to the probability of observing the last category in
#' state one.
#' }
#'
#' @param T The length of the observed sequence to be simulated for each
#'   subject. To only simulate subject specific transition probability matrices
#'   gamma and emission distributions (and no data), set \code{T} to 0.
#' @param n The number of subjects for which data is simulated.
#' @param m The number of hidden states in the HMM used for
#'   simulating data.
#' @param pr The number of categories of the simulated observations.
#' @param gamma A matrix with \code{m} rows and \code{m} columns containing the
#'   average population transition probability matrix used for simulating the
#'   data. That is, the probability to switch from hidden state \emph{i} (row
#'   \emph{i}) to hidden state \emph{j} (column  \emph{j}).
#' @param emiss_distr A matrix with \code{m} rows and \code{pr} columns
#'   containing the average population emission distribution of the
#'   (categorical) observations given the hidden states. That is, the
#'   probability of observing category \emph{k} (column \emph{k}) in state
#'   \emph{i} (row \emph{i}).
#' @param beta List of regressoin parameters to predict \code{gamma} and/or
#'   \code{emiss_distr} in combination with \code{xxS} using multinomial logisitc
#'   regression. The list is composed of two elements. The first element
#'   consists of a matrix with regression paremeters used to to predict the
#'   transition probability matrix \code{gamma}. The second element consist of a
#'   matrix with regression paremeters used to to predict the emission
#'   distribution \code{emiss_distr} of the dependent variable. Note that for
#'   both \code{gamma} and/or \code{emiss_distr}, the first element in each row
#'   is used as reference category in the multinomial logistic regression, and
#'   no regression parameters can be specified for these categories. Hence, the
#'   first element in the list \code{beta} to predict \code{gamma} consist of a
#'   matrix with the number of rows equal to \code{m} and the number of columns
#'   equal to \code{m} - 1. See \emph{details} for more information.  Note that
#'   if \code{beta} is specified, \code{xxS} has to be specified as well. If
#'   \code{beta} is omitted completely, \code{beta} defaults to NULL, resembling
#'   no prediction of \code{gamma} or \code{emiss_distr} using covariates.
#' @param xxS List of covariates to predict \code{gamma} and/or
#'   \code{emiss_distr} using the regression parameters specified in
#'   \code{beta}. The list is composed of two elements. The first element is
#'   used to predict the transition matrix. The second elements is used to
#'   predict the emission distribution of the dependent variable. Both elements
#'   are a matrix, with the number of rows equal to \code{n} the number of
#'   subjects to be simulated. If \code{xx} is omitted completely, \code{xx}
#'   defaults to NULL, resembling no covariates. Specific elements in the list
#'   can also be left empty (i.e., set to NULL) to signify that either the
#'   transition probability matrix or a specific emission distribution is not
#'   predicted by covariates.
#' @param var_gamma An integer denoting the variance between subjects in the
#'   transition probability matrix. Note that this value corresponds to the
#'   variance of the parameters of the multinomial distribution (i.e., the
#'   intercepts of the regression equation of the multinomial distribution used
#'   to sample the transition probability matrix), see details below. In
#'   addition, only one variance value can be specified for the complete
#'   transition probability matrix, hence the variance is assumed fixed across
#'   all components. The default equals 1, which corresponds to quite some
#'   variation between subjects. A less extreme value would be 0.5. If one wants
#'   to simulate data from exactly the same HMM for all subjects, var_gamma
#'   should be set to 0.
#' @param var_emiss An integer denoting the variance between subjects in the
#'   emission distribution. Note that this value corresponds to the variance of
#'   the parameters of the multinomial distribution (i.e., the intercepts of the
#'   regression equation of the multinomial distribution used to sample the
#'   components of the emission distribution), see details below.  In addition,
#'   only one variance value can be specified for the complete emission
#'   distribution, hence the variance is assumed fixed across all components.
#'   The default equals 1, which corresponds to quite some variation between
#'   subjects. A less extreme value would be 0.5. If one wants to simulate data
#'   from exactly the same HMM for all subjects, var_emmis should be set to
#'   0.
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
#'   is repeated over the rows (with the number of repeats equal to the lenght
#'   of the simulated hidden state sequence \code{T} for each subject).}
#'   \item{\code{obs}}{A matrix containing the simulated observed outputs, with
#'   one row per simulated observation per subject. The first column indicates
#'   subject id number. The second column contains the simulated observation
#'   sequence, consecutively for all subjects. Hence, the id number is repeated
#'   over rows (with the number of repeats equal to the lenght of the simulated
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
#'  @seealso \code{\link{mHMM_mnl}} for analyzing multilevel hidden Markov data.
#'
#'
#'   @examples
#' # simulating data for 10 subjects with each 100 observations
#' T <- 100
#' n <- 10
#' m <- 3
#' pr <- 4
#' gamma <- matrix(c(0.8, 0.1, 0.1,
#'                   0.2, 0.7, 0.1,
#'                   0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
#' emiss_distr <- matrix(c(0.5, 0.5, 0.0, 0.0,
#'                         0.1, 0.1, 0.8, 0.0,
#'                         0.0, 0.0, 0.1, 0.9), nrow = m, ncol = pr, byrow = TRUE)
#' data1 <- sim_mHMM(T = T, n = n, m = m, pr = pr, gamma = gamma, emiss_distr = emiss_distr,
#'                   var_gamma = 1, var_emiss = 1)
#'
#'
#' # simulating subject specific transition probability matrices and emission distributions only
#' T <- 0
#' n <- 5
#' m <- 3
#' pr <- 4
#' gamma <- matrix(c(0.8, 0.1, 0.1,
#'                   0.2, 0.7, 0.1,
#'                   0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
#' emiss_distr <- matrix(c(0.5, 0.5, 0.0, 0.0,
#'                         0.1, 0.1, 0.8, 0.0,
#'                         0.0, 0.0, 0.1, 0.9), nrow = m, ncol = pr, byrow = TRUE)
#' data2 <- sim_mHMM(T = T, n = n, m = m, pr = pr, gamma = gamma, emiss_distr = emiss_distr,
#'                   var_gamma = 1, var_emiss = 1)
#' data2
#'
#' data3 <- sim_mHMM(T = T, n = n, m = m, pr = pr, gamma = gamma, emiss_distr = emiss_distr,
#'                   var_gamma = .5, var_emiss = .5)
#' data3
#' @export

sim_mHMM <- function(T, n, m, pr, gamma, emiss_distr, beta = NULL, xxS = NULL,
                     var_gamma = 1, var_emiss = 1, return_ind_par = FALSE){
  if (dim(gamma)[1] != m){
    stop(paste("The transiton probability matrix gamma should be a", m, "by", m, "matrix."))
  }
  if (dim(gamma)[2] != m){
    stop(paste("The transiton probability matrix gamma should be a", m, "by", m, "matrix."))
  }
  if (dim(emiss_distr)[1] != m){
    stop(paste("The number of rows of the emission distribution matrix should be equal to the number of states, which is", m, "."))
  }
  if (dim(emiss_distr)[2] != pr){
    stop(paste("The number of columns of the emission distribution matrix should be equal to the number of observable categories, which is", pr, "."))
  }
  states <- matrix(ncol = 2, nrow = T*n)
  states[,1] <- rep(1:n, each = T)
  obs <- matrix(ncol = 2, nrow = T*n)
  obs[,1] <- rep(1:n, each = T)
  sub_gamma <- rep(list(NULL), n)
  sub_emiss <- rep(list(NULL), n)
  mnl_gamma <- prob_to_int(gamma)
  mnl_emiss <- prob_to_int(emiss_distr)
  for(j in 1:n){
    sub_gamma[[j]] <- int_to_prob(mnl_gamma +
                                    rnorm(n = m * (m-1), mean = 0, sd = sqrt(var_gamma)))
    sub_emiss[[j]] <- int_to_prob(mnl_emiss +
                                    rnorm(n = m * (pr-1), mean = 0, sd = sqrt(var_emiss)))
    if(T != 0){
      init <- solve(t(diag(m) - sub_gamma[[j]] + 1), rep(1, m))
      states[((j-1) * T + 1), 2] <- sample(x = 1:m, size = 1, prob = init)
      obs[((j-1) * T + 1), 2] <- sample(x = 1:pr, size = 1, prob = sub_emiss[[j]][states[((j-1) * T + 1), 2],])
      for(i in 2:T){
        states[((j-1) * T + i), 2] <- sample(x = 1:m, size = 1, prob = sub_gamma[[j]][states[((j-1) * T + i - 1), 2],])
        obs[((j-1) * T + i), 2] <- sample(x = 1:pr, size = 1, prob = sub_emiss[[j]][states[((j-1) * T + i), 2],])
      }
    }
  }
  if (return_ind_par == FALSE & T != 0){
    return(list(states = states, obs = obs))
  } else if (return_ind_par == TRUE & T != 0){
    return(list(states = states, obs = obs, subject_gamma = sub_gamma, subject_emmis = sub_emiss))
  } else if (T == 0){
    return(list(subject_gamma = sub_gamma, subject_emmis = sub_emiss))
  }
}
