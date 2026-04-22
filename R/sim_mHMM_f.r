#' Simulate data using a multilevel hidden Markov model
#'
#' \code{sim_mHMM_f} simulates data for multiple subjects, for which the data
#' have a continuous (i.e., normally distributed) observations that follow a
#' hidden Markov model (HMM) with a fixed over subjects emission distribution
#' and a multilevel structure imposed on the transition probabilities. The
#' multilevel structure implies that each subject is allowed to have its own set
#' of parameters, and that the parameters at the subject level (level 1) are
#' tied together by a population distribution at level 2 for each of the
#' corresponding parameters. The shape of the population distribution for each
#' of the parameters is a normal distribution. In addition to (natural and/or
#' unexplained) heterogeneity between subjects, the subjects parameters can also
#' depend on a covariate.
#'
#' In simulating the data, having a multilevel structure means that the
#' parameters for each subject are sampled from the population level
#' distribution of the corresponding parameter. The user specifies the
#' population distribution for each parameter: the average population transition
#' probability matrix and its variance. For now, the variance of the mean population
#' parameters is assumed fixed for all components of the transition probability
#' matrix.
#'
#' One can simulate multivariate data. That is, the hidden states depend on more
#' than 1 observed variable simultaneously. The distributions of multiple
#' dependent variables for multivariate data are assumed to be independent.
#'
#' Note that the subject specific initial state distributions (i.e., the
#' probability of each of the states at the first time point) needed to simulate
#' the data are obtained from the stationary distributions of the subject
#' specific transition probability matrices gamma.
#'
#' \code{beta}: As the first element in each row of \code{gamma} is used as
#' reference category in the Multinomial logistic regression, the first matrix
#' in the list \code{beta} used to predict transition probability matrix
#' \code{gamma} has a number of rows equal to \code{m} and the number of columns
#' equal to \code{m} - 1. The first element in the first row corresponds to the
#' probability of switching from state one to state two. The second element in
#' the first row corresponds to the probability of switching from state one to
#' state three, and so on. The last element in the first row corresponds to the
#' probability of switching from state one to the last state.
#'
#'
#' @inheritParams mHMM_f
#' @param n_t Numeric vector with length 1 denoting the length of the observed
#'   sequence to be simulated for each subject. To only simulate subject
#'   specific transition probability matrices gamma and emission distributions
#'   (and no data), set \code{t} to 0.
#' @param n Numeric vector with length 1 denoting the number of subjects for
#'   which data is simulated.
#' @param m The argument \code{m} is deprecated; please specify using the input
#'   parameter \code{gen}.
#' @param n_dep The argument \code{n_dep} is deprecated; please specify using
#'   the input parameter \code{gen}.
#' @param start_state Optional numeric vector with length 1 denoting in which
#'   state the simulated state sequence should start. If left unspecified, the
#'   simulated state for time point 1 is sampled from the initial state
#'   distribution (which is derived from the transition probability matrix
#'   gamma).
#' @param data_distr A character vector of length 1 denoting the distribution
#'   adopted for the data given the hidden states. In this version of the
#'   package, it can only take the default value 'continuous', standing for the
#'   continuous observations following a normal distribution.
#' @param gamma A matrix with \code{m} rows and \code{m} columns containing the
#'   average population transition probability matrix used for simulating the
#'   data. That is, the probability to switch from hidden state \emph{i} (row
#'   \emph{i}) to hidden state \emph{j} (column  \emph{j}).
#' @param emiss_distr A list with \code{n_dep} elements containing the fixed
#'   over subjects emission distribution(s) of the observations given the hidden
#'   states for each of the dependent variables. Each element is a matrix with
#'   \code{m} rows and 2 columns; the first column denoting the mean of state
#'   \emph{i} (row \emph{i}) and the second column denoting the standard
#'   deviation of state \emph{i} (row \emph{i}) of the Normal distribution.
#' @param xx_vec List of 1 vector containing the covariate to predict the
#'   transition probability matrix \code{gamma} using the regression parameter
#'   specified in \code{beta} (see below). At this point, it is only possible to
#'   use one covariate for both \code{gamma}. The number of observations in the
#'   vectors should be  equal to the number of subjects to be simulated
#'   \code{n}. If \code{xx_vec} is omitted completely, \code{xx_vec} defaults to
#'   NULL, resembling no covariates at all.
#' @param beta List of 1 matrices containing the regression parameter to predict
#'   \code{gamma} in combination with \code{xx_vec} using (Multinomial logistic)
#'   regression. One regression parameter is specified for each element in
#'   \code{gamma}, with the following exception: the first element in each row
#'   of \code{gamma} is used as reference category in the Multinomial logistic
#'   regression. As such, no regression parameters can be specified for these
#'   parameters. Hence, the first element in the list \code{beta} to predict
#'   \code{gamma} consist of a matrix with the number of rows equal to \code{m}
#'   and the number of columns equal to \code{m} - 1. See \emph{details} for
#'   more information.
#'
#'   Note that if \code{beta} is specified, \code{xx_vec} has to be specified as
#'   well. If \code{beta} is omitted completely, \code{beta} defaults to NULL,
#'   resembling no prediction of \code{gamma} using a covariate.
#'
#' @param var_gamma Either a numeric vector with length 1 or a matrix of
#'   (\code{m} by \code{m} - 1) elements denoting the amount of variance
#'   between subjects in the transition probability matrix. Note that the
#'   value(s) correspond to the variance of the parameters of the Multinomial
#'   distribution (i.e., the intercepts of the regression equation of the
#'   Multinomial distribution used to sample the transition probability
#'   matrix), see details below. Also note that if only one variance value is
#'   provided, it will be adopted for the complete transition probability
#'   matrix, hence the variance is assumed fixed across all components. The
#'   default equals 0.1, which corresponds to little variation between
#'   subjects. If one wants to simulate data from exactly the same HMM for all
#'   subjects, var_gamma should be set to 0. Note that if data for only 1
#'   subject is simulated (i.e., n = 1), \code{var_gamma} is set to 0.
#' @param return_ind_par A logical scalar. Should the subject specific
#'   transition probability matrix \code{gamma} be returned by the function
#'   (\code{return_ind_par = TRUE}) or not (\code{return_ind_par = FALSE}). The
#'   default equals \code{return_ind_par = FALSE}.

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
#' }
#'
#'
#' @seealso \code{\link{mHMM_f}} for analyzing multilevel hidden Markov data.
#'
#'
#' @examples
#'
#' ## Example on multivariate continuous data
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
#' emiss_distr <- list(matrix(c( 50, 10,
#'                               100, 10,
#'                               150, 10), nrow = m, byrow = TRUE),
#'                     matrix(c(5, 2,
#'                              10, 5,
#'                              20, 3), nrow = m, byrow = TRUE))
#'
#' data_cont <- sim_mHMM_f(n_t = n_t, n = n, data_distr = 'continuous',
#'                       gen = list(m = m, n_dep = n_dep),
#'                       gamma = gamma, emiss_distr = emiss_distr,
#'                       var_gamma = .5)
#'
#' head(data_cont$states)
#' head(data_cont$obs)
#'
#'
#'



#' @export

sim_mHMM_f <- function(n_t, n, data_distr = 'continuous', gen, gamma, emiss_distr, start_state = NULL,
                     xx_vec = NULL, beta = NULL, var_gamma = 0.1, return_ind_par = FALSE, m, n_dep){
  if(data_distr != 'continuous'){
    stop("this developers version can only be used to simulate observations that have a Normal (i.e., Gaussian) emission distribution")
  }
  if(!missing(m)){
    warning("The argument m is deprecated; please specify using the input parameter gen.")
  }
  if(!missing(n_dep)){
    warning("The argument n_dep is deprecated; please specify using the input parameter gen.")
  }

  if(!missing(gen)){
    if(sum(objects(gen) %in% "m") != 1 | sum(objects(gen) %in% "n_dep") != 1 ){
      stop("The input argument gen should contain the elements m and n_dep.")
    }
    m <- gen$m
    n_dep <- gen$n_dep
  }
  if(m == 1){
    warning("If the number of states m is set to 1, variance between subjects in the transition probability matrix is ignored, and cannot be explained by covariates.")
  }
  if(missing(m) & missing(gen)){
    stop("Please specify the number of hidden states m via the input parameter gen.")
  }
  if(missing(n_dep) & missing(gen)){
    n_dep <- 1
    warning("Please specify the number of dependent variables n_dep via the input parameter gen. Model now assums n_dep = 1")
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
  for(q in 1:n_dep){
    if (dim(emiss_distr[[q]])[1] != m){
      stop(paste("The number of rows of emission distribution matrix in element", q, "should be
             equal to the number of states, which is", m, "."))
    }
    if(data_distr == 'continuous'){
      if (dim(emiss_distr[[q]])[2] != 2){
        stop(paste("For continuous data, the number of columns of the emission distribution matrix should be 2, where the first column denotes the state dependent mean and the
                  second column the state dependent standard deviation of the Normal emission distribution. See emission distribution in element", q, "."))
      }
    }
  }
  if((is.null(xx_vec) & !is.null(beta)) | (!is.null(xx_vec) & is.null(beta))){
    stop("Either only xx_vec or only beta is specified. Please specify both 1) the values for the covariate
         in xx_vec and 2) the values of the regression parameters in beta, to allow correct simulation of the
         data.")
  }
  if(!is.null(xx_vec)){
    if(!is.list(xx_vec) |( length(xx_vec) > 1)){
      stop("The format of xx_vec should be a list with - for this developers version - only one element")
    }
  }
  if(!is.null(beta)){
    if(!is.list(beta) |( length(beta) > 1)){
      stop("The format of beta should be a list with - for this developers version - only one element")
    }
  }
  if(!is.null(xx_vec)){
    # extend to all 1 + n_dep
    if((!is.null(xx_vec[[1]]) & length(xx_vec[[1]]) != n)){
      stop("The length of the vectors in xx_vec should be equal to the number of subjects to be simulated,
         set in n, if xx_vec is not set to NULL.")
    }
  }
  if (!is.null(beta)){
    if (!is.null(beta[[1]])){
      if ((dim(beta[[1]])[1] != (m)) | (dim(beta[[1]])[2] != (m-1))){
        stop(paste("The first element of beta to predict the transiton probability matrix gamma should be a m (", m, " ) by m - 1 (", m - 1, ") matrix."))
      }
    }
  }
  if(is.null(xx_vec)){
    xx_vec <- rep(list(NULL), 1)
    xx_vec[[1]] <- rep(1,n)
  }
  if(is.null(beta)){
    beta <- rep(list(NULL), 1)
    beta[[1]] <- matrix(0, ncol = m - 1, nrow = m)
  }

  # If only 1 subject
  if(n == 1){
    var_gamma <- 0
    var_emiss <- rep(0, n_dep)
  }

  # If a single value of var_gamma specified, use for all transitions
  if(length(var_gamma) == 1){
    var_gamma <- matrix(rep(var_gamma, m*(m-1)),nrow = m, byrow = TRUE)
    # warning("A single value of var_gamma was provided, which will be used for all states.")
  } else if(is.matrix(var_gamma)){
    if (dim(var_gamma)[1] != m){
      stop(paste("The between-subject variance matrix for the transition distribution should be a", m, "by", m-1, "matrix."))
    }
    if (dim(var_gamma)[2] != m-1){
      stop(paste("The between-subject variance matrix for the transition distribution should be a", m, "by", m-1, "matrix."))
    }
  }


  #############
  # Simulating the data ---------------------
  #############

  states <- matrix(ncol = 2, nrow = n_t*n)
  states[,1] <- rep(1:n, each = n_t)
  obs <- matrix(ncol = 1 + n_dep, nrow = n_t*n)
  obs[,1] <- rep(1:n, each = n_t)
  if(m > 1){
    sub_gamma <- rep(list(NULL), n)
    mnl_gamma <- prob_to_int(gamma)
  } else if (m == 1){
    sub_gamma <- rep(list(matrix(1)), n)
  }


  for(j in 1:n){
    if(m > 1){
      sub_gamma[[j]] <- int_to_prob(mnl_gamma + xx_vec[[1]][j] * beta[[1]] +
                                      rnorm(n = m * (m-1), mean = 0, sd = sqrt(as.numeric(var_gamma))))
    }

    if(n_t != 0){
      init <- solve(t(diag(m) - sub_gamma[[j]] + 1), rep(1, m))
      if (is.null(start_state)){
        states[((j-1) * n_t + 1), 2] <- sample(x = 1:m, size = 1, prob = init)
      } else {
        states[((j-1) * n_t + 1), 2] <- start_state
      }
      for(i in 1:n_dep){
        obs[((j-1) * n_t + 1), (1+i)] <- rnorm(1, mean = emiss_distr[[i]][states[((j-1) * n_t + 1), 2],1],
                                                 sd = emiss_distr[[i]][states[((j-1) * n_t + 1), 2],2])
      }
      for(t in 2:n_t){
        states[((j-1) * n_t + t), 2] <- sample(x = 1:m, size = 1, prob = sub_gamma[[j]][states[((j-1) * n_t + t - 1), 2],])
        for(i in 1:n_dep){
          obs[((j-1) * n_t + t), (1+i)] <- rnorm(1, mean = emiss_distr[[i]][states[((j-1) * n_t + t), 2],1],
                                                 sd = emiss_distr[[i]][states[((j-1) * n_t + t), 2],2])
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
    return(list(states = states, obs = obs, subject_gamma = sub_gamma))
  } else if (n_t == 0){
    return(list(subject_gamma = sub_gamma))
  }
}
