#' Multilevel hidden markov model using Bayesian estimation
#'
#' \code{mHMM_mnl} analyses (intense longitudinal) data from multiple subjects
#' using a multilevel hidden Markov model. By using a multilevel framework, one
#' general 'population' HMM is estimated, while heterogeneity between subjects
#' is accommodated. The function can handle covariates at the subject level
#' (unlimited number), uses a hybrid metropolis within gibs sampler, and
#' performs the forward backward algorithm for all subjects in a sequential
#' manner. Can handle varying observation length over subjects.
#'
#' @param s_data A matrix containing the observations to be modelled, with one
#'   row per observation. In the \code{s_data} matrix, the first column
#'   indicates subject id number. Hence, the id number is repeated over rows
#'   equal to the number of observations for that subject. The number of
#'   observations can vary over subjects. The subsequent columns contain the
#'   dependent variable(s). The total number of rows are equal to the sum over
#'   the number of observations of each subject, and the number of columns are
#'   equal to the number of dependent variables (\code{n_dep}) + 1.
#' @param gen List containing the following elements:
#'   \itemize{\item{\code{m}: numeric vector with length 1 denoting the number
#'   of states}
#'   \item{\code{n_dep}: numeric vector with length 1 denoting the
#'   number of dependent variables}
#'   \item{\code{q_emiss}: numeric vector with length \code{n_dep} denoting the
#'   number of observed categories for the categorical  emission distribution of
#'   each dependent variable.}}
#' @param xx List of covariates. Number of elements in the list is equal to 1 +
#'   \code{n_dep} (i.e., the number of dependent variables). The first element
#'   in the list is used to predict the transition matrix. Subsequent elements
#'   in the list are used to predict the emission distribution of (each of) the
#'   dependent variable(s). Each element in the list is a matrix, with the
#'   number of rows equal to the number of subjects. The first column \emph{has
#'   to}  represent the intercept, that is, a column only consisting of ones.
#'   Subsequent columns correspond to covariates used to predict the transition
#'   matrix / emission distribution. Covariates can either be dichotomous or
#'   continuous variables. Dichotomous variables have to be coded as 0/1
#'   variables. Categroical or factor variables can as yet not be used as
#'   predictor covariates. The user can however break up the categorical
#'   variable in multiple dummy variables (i.e., dichotomous variables), which
#'   can be used simultaneously in the analysis. Continuous predictors are
#'   automatically centered. That is, the mean value of the covariate is
#'   subtracted from all values of the covariate such that the new mean equals
#'   zero. This is done such that the presented probabilities in the output
#'   (i.e., for the population transiton probability matrix and population
#'   emission probabilities) correspond to the predicted probabilities at the
#'   average value of the covariate(s). MOVE PART TO DETAILS
#'
#'   If \code{xx} is omitted completely, \code{xx} defaults to NULL, resembling
#'   no covariates. Specific elements in the list can also be left empty (i.e.,
#'   set to NULL) to signify that either the transition probability matrix or a
#'   specific emission distribution is not predicted by covariates.
#' @param start_val List containing the start values for the transition probability
#'   matrix gamma and the emission distribuition(s). The first element of the list
#'   contains a matrix with the start values for gamma. The subsequent elements
#'   contain matrices for the start values for the conditional distribution(s):
#'   one element containing one matrix for each emission distribtution of a
#'   dependent variable.
#' @param gamma_sampler List containing start values for mle estimates of pooled
#'   data for gamma: \code{gamma_int_mle0}, \code{gamma_scalar}, and weight for the overall ll in the
#'   fractional likelihood, \code{gamma_w}.
#' @param emiss_sampler List containing start values for mle estimates of pooled
#'   data for emiss, emiss_int_mle0, emiss_scalar and weight for the overall ll in the
#'   fractional likelihood, emiss_w
#' @param gamma_hyp_prior List containing \code{gamma_mu0}, \code{gamma_K0}, \code{gamma_nu} and \code{gamma_V}, see details below
#' @param emiss_hyp_prior List containing \code{emiss_mu0}, \code{emiss_K0}, \code{emiss_nu} and \code{emiss_V}, see details below
#' @param mcmc List of MCMC argument, containing the following elements:
#'   \code{J} number of iterations of the MCMC algorithm. \code{burn_in} Burn in
#'   period for the MCMC algorithm
#' @param return_path A logical scalar. Should the sampled state sequence obtained at
#'   each iteration and for each subject be returned by the function
#'   (\code{sample_path = TRUE}) or not (\code{sample_path = FALSE}). This is
#'   quite a large object! Can be used for local decoding purposes.
#'
# @details Here are the details of the function
# \subsection{Equations}{
#   The following hyper prior is used on each row \eqn{i} of gamma (\eqn{\Gamma_i}) at the population level:
#   \deqn{\Gamma_i ~ (idd) MNL(gamma_int_i)}
#   \deqn{gamma_int_i ~ N(gamma_mu_int_bar, gamma_V_int)}
#   \deqn{gamma_mu_int_bar ~ N(gamma_mu0, \frac{1}{gamma_K0} * gamma_V_int)}
#   \deqn{gamma_V_int ~ IW(gamma_nu, gamma_V)}
#   where \eqn{\alpha_i} denotes \code{gamma_int_i}, \eqn{\bar{\mu}_\alpha} denotes \code{gamma_mu_int_bar},  \code{} \code{} \code{} \code{} \code{}
#
#   The following hyper prior is used for the conditional emission distribution of each state at the population level:
#   \deqn{\theta_i ~ (idd) MNL(emiss_int_i)}
#   \deqn{emiss_int ~ N(emiss_mu_int_bar, emiss_V_int)}
#   \deqn{emiss_mu_int_bar ~ N(emiss_mu0, (1 / emiss_K0) * emiss_V_int)}
#   \deqn{emiss_V_int 					~ IW(emiss_nu, emiss_V)}
#
# check thesis for correct notation and stuff
# }
#
# \subsection{Covariates}{
#   Here is some info on using covariates ... (to be completed)
# }
#
#' @return \code{mHMM_mnl} returns an object (list) of the class \code{mHMM}.
#'   The list contains the following components:
#'   \describe{
#'   \item{\code{PD_subj}}{A list with one matrix per subject containing the xx,
#'   log likelihood}
#'   \item{\code{emiss_prob_bar}}{A matrix containing XX, with
#'   one row per XX. The columns ... }
#'   \item{\code{emiss_int_bar}}{A matrix containing XX, with
#'   one row per XX. The columns ... }
#'   \item{\code{emiss_int_subj}}{A matrix containing XX, with
#'   one row per XX. The columns ... }
#'   \item{\code{emiss_naccept}}{A matrix containing XX, with
#'   one row per XX. The columns ... }
#'   \item{\code{gamma_prob_bar}}{A matrix containing XX, with
#'   one row per XX. The columns ... }
#'   \item{\code{gamma_int_bar}}{A matrix containing XX, with
#'   one row per XX. The columns ... }
#'   \item{\code{gamma_int_subj}}{A matrix containing XX, with
#'   one row per XX. The columns ... }
#'   \item{\code{gamma_naccept}}{A matrix containing XX, with
#'   one row per XX. The columns ... }
#'   \item{\code{input}}{Overview of used input specifications. Namely, ...}
#'   \item{\code{sample_path}}{Only returned if .. }
#' }
#'
#' @seealso \code{\link{sim_mHMM}} for simulating multilevel hidden Markov data
#'   and \code{\link{vit_mHMM}} for obtaining the most likely hidden state
#'   sequence for each subject using the Viterbi algorithm.
#'
#' @examples
#' ###### Example on package data
#' # specifying general model properties:
#' m <- 2
#' n_dep <- 4
#' q_emiss <- c(3, 2, 3, 2)
#'
#' # specifying starting values
#' start.TM <- diag(.8, m)
#' start.TM[lower.tri(start.TM) | upper.tri(start.TM)] <- .2
#' start.EM <- list(matrix(c(0.9, 0.05, 0.05, 0.05, 0.05, 0.9), byrow = TRUE,
#'                         nrow = m, ncol = q_emiss[1]), # vocalizing patient
#'                  matrix(c(0.9, 0.1, 0.9, 0.1), byrow = TRUE, nrow = m,
#'                         ncol = q_emiss[2]), # looking patient
#'                  matrix(c(0.05, 0.05, 0.9, 0.9, 0.05, 0.05), byrow = TRUE,
#'                         nrow = m, ncol = q_emiss[3]), # vocalizing therapist
#'                  matrix(c(0.9, 0.1, 0.9, 0.1), byrow = TRUE, nrow = m,
#'                         ncol = q_emiss[4])) # looking therapist
#'
#' # Run a model without covariate(s):
#' out1 <- mHMM_mnl(s_data = nonverbal, gen = list(m = m, n_dep = n_dep,
#'                 q_emiss = q_emiss), start_val = c(list(start.TM), start.EM),
#'                 mcmc = list(J = 11, burn_in = 5))
#'
#' # plot the posterior densities for the transition and emission probabilities
#' plot(out1, component = "gamma")
#' plot(out1, component = "emiss", dep = 1)
#' plot(out1, component = "emiss", dep = 2)
#' plot(out1, component = "emiss", dep = 3)
#' plot(out1, component = "emiss", dep = 4)
#'
#' # Run a model including a covariate. Here, the covariate (standardized CDI
#' # change) predicts the emission distribution for each of the 4 dependent
#' # variables:
#' n_subj <- 10
#' xx <- rep(list(matrix(1, ncol = 1, nrow = n_subj)), (n_dep + 1))
#' for(i in 2:(n_dep + 1)){
#'   xx[[i]] <- cbind(xx[[i]], nonverbal_cov$std_CDI_change)
#' }
#' out2 <- mHMM_mnl(s_data = nonverbal, xx = xx, gen = list(m = m, n_dep = n_dep,
#'                 q_emiss = q_emiss), start_val = c(list(start.TM), start.EM),
#'                 mcmc = list(J = 11, burn_in = 5))
#'
#'
#' ###### Example on simulated data
#' # Simulate data for 10 subjects with each 100 observations:
#' n_t <- 100
#' n <- 10
#' m <- 2
#' q_emiss <- 3
#' gamma <- matrix(c(0.8, 0.2,
#'                   0.3, 0.7), ncol = m, byrow = TRUE)
#' emiss_distr <- matrix(c(0.5, 0.5, 0.0,
#'                         0.1, 0.1, 0.8), nrow = m, ncol = q_emiss, byrow = TRUE)
#' data1 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
#'                   emiss_distr = emiss_distr, var_gamma = .5, var_emiss = .5)
#'
#' # Specify remaining required anaylsis input (for now, use simulation input as
#' # starting values):
#' n_dep <- 1
#' q_emiss <- 3
#'
#' # Run the model on the simulated data:
#' out3 <- mHMM_mnl(s_data = data1$obs, gen = list(m = m, n_dep = n_dep,
#'                 q_emiss = q_emiss), start_val = c(as.vector(t(emiss_distr)),
#'                 as.vector(t(gamma))), mcmc = list(J = 11, burn_in = 5))
#'
#' @export
#'
#'

mHMM_mnl <- function(s_data, gen, xx = NULL, start_val, gamma_sampler = NULL, emiss_sampler = NULL,
                     gamma_hyp_prior = NULL, emiss_hyp_prior = NULL, mcmc, return_path = FALSE){

  # Initialize data -----------------------------------
  # dependent variable(s), sample size, dimensions gamma and conditional distribuiton
  n_dep			 <- gen$n_dep
  dep_labels <- colnames(s_data[,2:(n_dep+1)])
  id         <- unique(s_data[,1])
  n_subj     <- length(id)
  subj_data  <- rep(list(NULL), n_subj)
  for(s in 1:n_subj){
    subj_data[[s]]$y <- as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep)
  }
  ypooled    <- n <- NULL
  n_vary     <- numeric(n_subj)
  m          <- gen$m
  q_emiss 		 <- gen$q_emiss
  emiss_int_mle <- rep(list(NULL), n_dep)
  emiss_mhess   <- rep(list(NULL), n_dep)
  for(q in 1:n_dep){
    emiss_int_mle[[q]] <- matrix(, m, (q_emiss[q] - 1))
    emiss_mhess[[q]] <- matrix(, (q_emiss[q] - 1) * m, (q_emiss[q] - 1))
  }
  for(s in 1:n_subj){
    ypooled   <- rbind(ypooled, subj_data[[s]]$y)
    n         <- dim(subj_data[[s]]$y)[1]
    n_vary[s] <- n
    subj_data[[s]]	<- c(subj_data[[s]], n = n, list(gamma_converge = numeric(m), gamma_int_mle = matrix(, m, (m - 1)),
                                                    gamma_mhess = matrix(, (m - 1) * m, (m - 1)), emiss_converge =
                                                      rep(list(numeric(m)), n_dep), emiss_int_mle = emiss_int_mle, emiss_mhess = emiss_mhess))
  }
  n_total 		<- dim(ypooled)[1]

  # covariates
  n_dep1 <- 1 + n_dep
  nx <- numeric(n_dep1)
  if (is.null(xx)){
    xx <- rep(list(matrix(1, ncol = 1, nrow = n_subj)), n_dep1)
    nx[] <- 1
  } else {
    for(i in 1:n_dep1){
      if (is.null(xx[[i]])){
        xx[[i]] <- matrix(1, ncol = 1, nrow = n_subj)
        nx[i] <- 1
      } else {
        nx[i] <- ncol(xx[[i]])
        if (sum(xx[[i]][,1] != 1)){
          stop("If xx is specified, the first column in each element of xx has to represent the intercept. That is, a column that only consists of the value 1")
        }
        if(nx[i] > 1){
          for(j in 2:nx[i]){
            if(is.factor(xx[[i]][,j])){
              stop("Factors currently cannot be used as covariates, see help file for alternatives")
            }
            if((length(unique(xx[[i]][,j])) == 2) & (sum(xx[[i]][,j] != 0 & xx[[i]][,j] !=1) > 0)){
              stop("Dichotomous covariates in xx need to be coded as 0 / 1 variables. That is, only conisting of the values 0 and 1")
            }
            if(length(unique(xx[[i]][,j])) > 2){
              xx[[i]][,j] <- xx[[i]][,j] - mean(xx[[i]][,j])
            }
          }
        }
      }
    }
  }

  # Initialize mcmc argumetns
  J 				<- mcmc$J
  burn_in			<- mcmc$burn_in


  # Initalize priors and hyper priors --------------------------------
  # Initialize gamma sampler
  if(is.null(gamma_sampler)) {
    gamma_int_mle0  <- gamma_mu_prop <- rep(0, m - 1)
    gamma_scalar    <- 2.93 / sqrt(m - 1)
    gamma_w         <- .1
  } else {
    gamma_int_mle0  <- gamma_sampler$gamma_int_mle0
    gamma_mu_prop   <- gamma_sampler$gamma_mu_prop
    gamma_scalar    <- gamma_sampler$gamma_scalar
    gamma_w         <- gamma_sampler$gamma_w
  }

  # Initialize emiss sampler
  if(is.null(emiss_sampler)){
    emiss_int_mle0 <- rep(list(NULL), n_dep)
    emiss_scalar	<- rep(list(NULL), n_dep)
    emiss_mu_prop 	<-  rep(list(NULL), n_dep)
    for(q in 1:n_dep){
      emiss_int_mle0[[q]]	<- emiss_mu_prop[[q]] <- rep(0, q_emiss[q] - 1)
      emiss_scalar[[q]] 	<- 2.93 / sqrt(q_emiss[q] - 1)
    }
    emiss_w		<- .1
  } else {
    emiss_int_mle0	<- emiss_sampler$emiss_int_mle0
    emiss_mu_prop 	<- emiss_sampler$emiss_mu_prop
    emiss_scalar 	<- emiss_sampler$emiss_scalar
    emiss_w    		<- emiss_sampler$emiss_w
  }

  # Initialize Gamma hyper prior
  if(is.null(gamma_hyp_prior)){
    gamma_mu0			<- matrix(0,nrow = nx, ncol = m - 1)
    gamma_K0			<- diag(1, nx[1])
    gamma_nu			<- 3 + m - 1
    gamma_V			  <- gamma_nu * diag(m - 1)
  } else {
    ###### BUILD in a warning / check if gamma_mu0 is a matrix when given, with  nrows equal to the number of covariates
    gamma_mu0			<- gamma_hyp_prior$gamma_mu0
    gamma_K0			<- gamma_hyp_prior$gamma_K0
    gamma_nu			<- gamma_hyp_prior$gamma_nu
    gamma_V			  <- gamma_hyp_prior$gamma_V
  }

  # Initialize Pr hyper prior
  # emiss_mu0: for each dependent variable, emiss_mu0 is a list, with one element for each state.
  # Each element is a matrix, with number of rows equal to the number of covariates (with the intercept being one cov),
  # and the number of columns equal to q_emiss[q] - 1.
  emiss_mu0	  <- rep(list(vector("list", m)), n_dep)
  emiss_nu	    <- rep(list(NULL), n_dep)
  emiss_V	    <- rep(list(NULL), n_dep)
  emiss_K0     <- rep(list(NULL), n_dep)
  if(is.null(emiss_hyp_prior)){
    for(q in 1:n_dep){
      for(i in 1:m){
        emiss_mu0[[q]][[i]]		<- matrix(0, ncol = q_emiss[q] - 1, nrow = nx[1 + q])
      }
      emiss_nu[[q]]		<- 3 + q_emiss[q] - 1
      emiss_V[[q]]		  <- emiss_nu[[q]] * diag(q_emiss[q] - 1)
      emiss_K0[[q]]		<- diag(1, nx[1 + q])
    }
  } else {
    for(q in 1:n_dep){
      # emiss_hyp_prior[[q]]$emiss_mu0 has to contain a list with lenght equal to m, and each list contains matrix with number of rows equal to number of covariates for that dep. var.
      # stil build in a CHECK, with warning / stop / switch to default prior
      emiss_mu0[[q]]	 <- emiss_hyp_prior[[q]]$emiss_mu0
      emiss_nu[[q]]	 <- emiss_hyp_prior[[q]]$emiss_nu
      emiss_V[[q]]		 <- emiss_hyp_prior[[q]]$emiss_V
      emiss_K0[[q]]	 <- diag(emiss_hyp_prior$emiss_K0, nx[1 + q])
    }
  }


  # Define objects used to store data in mcmc algorithm, not returned ----------------------------
  # overall
  c <- llk <- numeric(1)
  sample_path <- lapply(n_vary, dif_matrix, cols = J)
  trans <- rep(list(vector("list", m)), n_subj)

  # gamma
  gamma_mle_pooled <-gamma_int_mle_pooled <- gamma_mhess_pooled <- gamma_pooled_ll <- vector("list", m)
  gamma_c_int <- rep(list(matrix(, n_subj, (m-1))), m)
  gamma_mu_int_bar <- gamma_V_int <- vector("list", m)
  gamma_mu_prob_bar <- rep(list(numeric(m)), m)
  gamma_naccept <- matrix(0, n_subj, m)

  # emiss
  cond_y <- lapply(rep(n_dep, n_subj), nested_list, m = m)
  emiss_mle_pooled <- emiss_int_mle_pooled <- emiss_mhess_pooled <- emiss_pooled_ll <- rep(list(vector("list", n_dep)), m)
  emiss_c_int <- rep(list(lapply(q_emiss - 1, dif_matrix, rows = n_subj)), m)
  emiss_mu_int_bar <- emiss_V_int <- rep(list(vector("list", n_dep)), m)
  emiss_mu_prob_bar <- rep(list(lapply(q_emiss, dif_vector)), m)
  emiss_naccept <- rep(list(matrix(0, n_subj, m)), n_dep)


  # Define objects that are returned from mcmc algorithm ----------------------------
  # Define object for subject specific posterior density, put start values on first row
  PD 					  <- matrix(, nrow = J, ncol = sum(m * q_emiss) + m * m + 1)
  PD_emiss_names   <- paste("q", 1, "_emiss", rep(1:q_emiss[1], m), "_S", rep(1:m, each = q_emiss[1]), sep = "")
  if(n_dep > 1){
    for(q in 2:n_dep){
      PD_emiss_names <- c(PD_emiss_names, paste("q", q, "_emiss", rep(1:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q]), sep = ""))
    }
  }
  colnames(PD) 	<- c(PD_emiss_names, paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = ""), "LL")

  PD[1, ((sum(m * q_emiss) + 1)) :((sum(m * q_emiss) + m * m))] <- unlist(sapply(start_val, t))[1:(m*m)]
  PD[1, 1:((sum(m * q_emiss)))] <- unlist(sapply(start_val, t))[(m*m + 1): (m*m + sum(m * q_emiss))]

  PD_subj				<- rep(list(PD), n_subj)

  # Define object for population posterior density (probabilities and regression coefficients parameterization )
  gamma_prob_bar		<- matrix(, nrow = J, ncol = (m * m))
  colnames(gamma_prob_bar) <- paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = "")
  emiss_prob_bar			<- lapply(q_emiss * m, dif_matrix, rows = J)
  names(emiss_prob_bar) <- dep_labels
  for(q in 1:n_dep){
    colnames(emiss_prob_bar[[q]]) <- paste("Emiss", rep(1:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q]), sep = "")
  }
  gamma_int_bar				<- matrix(, nrow = J, ncol = ((m-1) * m))
  colnames(gamma_int_bar) <- paste("int_S", rep(1:m, each = m-1), "toS", rep(2:m, m), sep = "")
  if(nx[1] > 1){
    gamma_cov_bar				<- matrix(, nrow = J, ncol = ((m-1) * m) * (nx[1] - 1))
    colnames(gamma_cov_bar) <- paste( paste("cov", 1 : (nx[1] - 1), "_", sep = ""), "S", rep(1:m, each = (m-1) * (nx[1] - 1)), "toS", rep(2:m, m * (nx[1] - 1)), sep = "")
  } else{
    gamma_cov_bar <- "No covariates where used to predict the transition probability matrix"
  }
  emiss_int_bar			<- lapply((q_emiss-1) * m, dif_matrix, rows = J)
  names(emiss_int_bar) <- dep_labels
  for(q in 1:n_dep){
    colnames(emiss_int_bar[[q]]) <-  paste("int_Emiss", rep(2:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q] - 1), sep = "")
  }
  if(sum(nx[-1]) > n_dep){
    emiss_cov_bar			<- lapply((q_emiss-1) * m * (nx[-1] - 1 ), dif_matrix, rows = J)
    names(emiss_cov_bar) <- dep_labels
    for(q in 1:n_dep){
      if(nx[1 + q] > 1){
        colnames(emiss_cov_bar[[q]]) <-  paste( paste("cov", 1 : (nx[1 + q] - 1), "_", sep = ""), "emiss", rep(2:q_emiss[q], m * (nx[1 + q] - 1)), "_S", rep(1:m, each = (q_emiss[q] - 1) * (nx[1 + q] - 1)), sep = "")
      } else {
        emiss_cov_bar[[q]] <- "No covariates where used to predict the emission probabilities for this outcome"
      }
    }
  } else{
    emiss_cov_bar <- "No covariates where used to predict the emission probabilities"
  }

  # Define object for subject specific posterior density (regression coefficients parameterization )
  gamma_int_subj			<- rep(list(gamma_int_bar), n_subj)
  emiss_int_subj			<- rep(list(emiss_int_bar), n_subj)

  # Put starting values in place for fist run forward algorithm
  emiss_sep 			<- vector("list", n_dep)
  for(q in 1:n_dep){
    start <- c(0, q_emiss * m)
    emiss_sep[[q]] <- matrix(PD[1,(sum(start[1:q]) + 1):(sum(start[1:q]) + (m * q_emiss[q]))], byrow = TRUE, ncol = q_emiss[q], nrow = m)
  }
  emiss				  <- rep(list(emiss_sep), n_subj)
  gamma 			<- rep(list(matrix(PD[1,(sum(m*q_emiss) + 1):(sum(m * q_emiss) + m * m)], byrow = TRUE, ncol = m)), n_subj)
  delta 			<- rep(list(solve(t(diag(m) - gamma[[1]] + 1), rep(1, m))), n_subj)


  # Start analysis --------------------------------------------
  # Run the MCMC algorithm
  itime <- proc.time()[3]
  for (iter in 2 : J){

    # For each subject, obtain sampled state sequence with subject individual parameters ----------
    for(s in 1:n_subj){
      # Run forward algorithm, obtain subject specific forward proababilities and log likelihood
      forward				<- cat_Mult_HMM_fw(x = subj_data[[s]]$y, m = m, emiss = emiss[[s]], gamma = gamma[[s]], n_dep = n_dep, delta=NULL)
      alpha         <- forward$forward_p
      c             <- max(forward$la[, subj_data[[s]]$n])
      llk           <- c + log(sum(exp(forward$la[, subj_data[[s]]$n] - c)))
      PD_subj[[s]][iter, sum(m * q_emiss) + m * m + 1] <- llk

      # Using the forward probabilites, sample the state sequence in a backward manner.
      # In addition, saves state transitions in trans, and conditional observations within states in cond_y
      trans[[s]]					                  <- vector("list", m)
      sample_path[[s]][n_vary[[s]], iter] 	<- sample(1:m, 1, prob = c(alpha[, n_vary[[s]]]))
      for(t in (subj_data[[s]]$n - 1):1){
        sample_path[[s]][t,iter] 	              <- sample(1:m, 1, prob = (alpha[, t] * gamma[[s]][,sample_path[[s]][t + 1, iter]]))
        trans[[s]][[sample_path[[s]][t,iter]]]	<- c(trans[[s]][[sample_path[[s]][t, iter]]], sample_path[[s]][t + 1, iter])
      }
      for (i in 1:m){
        trans[[s]][[i]] <- c(trans[[s]][[i]], 1:m)
        trans[[s]][[i]] <- rev(trans[[s]][[i]])
        for(q in 1:n_dep){
          cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q], 1:q_emiss[q])
        }
      }
    }

    # The remainder of the mcmc algorithm is state specific
    for(i in 1:m){

# Obtain MLE of the covariance matrices and log likelihood of gamma and emiss at subject and population level -----------------
      # used to scale the propasal distribution of the RW Metropolis sampler

      # population level, transition matrix
      trans_pooled			  <- factor(c(unlist(sapply(trans, "[[", i)), c(1:m)))
      gamma_mle_pooled[[i]] 		<- optim(gamma_int_mle0, llmnl_int, Obs = trans_pooled, n_cat = m, method = "BFGS", hessian = TRUE, control = list(fnscale = -1))
      gamma_int_mle_pooled[[i]] <- gamma_mle_pooled[[i]]$par
      gamma_pooled_ll[[i]]			<- gamma_mle_pooled[[i]]$value

      # population level, conditional probabilities, seperate for each dependent variable
      for(q in 1:n_dep){
        cond_y_pooled					      <- numeric()
        ### MOET OOK ECHT BETER KUNNEN, eerst # cond_y_pooled				<- unlist(sapply(cond_y, "[[", m))
        for(s in 1:n_subj){
          cond_y_pooled             <- c(cond_y_pooled, cond_y[[s]][[i]][[q]])
        }
        emiss_mle_pooled[[i]][[q]]		  <- optim(emiss_int_mle0[[q]], llmnl_int, Obs = c(cond_y_pooled, c(1:q_emiss[q])), n_cat = q_emiss[q],
                                          method = "BFGS", hessian = TRUE, control = list(fnscale = -1))
        emiss_int_mle_pooled[[i]][[q]]	<- emiss_mle_pooled[[i]][[q]]$par
        emiss_pooled_ll[[i]][[q]]			<- emiss_mle_pooled[[i]][[q]]$value
      }

      # subject level, transition matrix
      for (s in 1:n_subj){
        wgt 				<- subj_data[[s]]$n / n_total
        gamma_out					<- optim(gamma_int_mle_pooled[[i]], llmnl_int_frac, Obs = c(trans[[s]][[i]], c(1:m)), n_cat = m,
                         pooled_likel = gamma_pooled_ll[[i]], w = gamma_w, wgt = wgt, method="BFGS", hessian = TRUE, control = list(fnscale = -1))
        if(gamma_out$convergence == 0){
          subj_data[[s]]$gamma_converge[i]									<- 1
          subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ]	<- mnlHess_int(int = gamma_out$par, Obs = c(trans[[s]][[i]], c(1:m)), n_cat =  m)
          subj_data[[s]]$gamma_int_mle[i,]									<- gamma_out$par
        } else {
          subj_data[[s]]$gamma_converge[i]									<- 0
          subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)),]	<- diag(m-1)
          subj_data[[s]]$gamma_int_mle[i,]									<- rep(0, m - 1)
        }
        # if this is first iteration, use MLE for current values RW metropolis sampler
        if (iter == 2){
          gamma_c_int[[i]][s,]		<- gamma_out$par
        }

        # subject level, conditional probabilities, seperate for each dependent variable
        for(q in 1:n_dep){
          emiss_out				<- optim(emiss_int_mle_pooled[[i]][[q]], llmnl_int_frac, Obs = c(cond_y[[s]][[i]][[q]], c(1:q_emiss[q])),
                             n_cat = q_emiss[q], pooled_likel = emiss_pooled_ll[[i]][[q]], w = emiss_w, wgt = wgt, method = "BFGS", hessian = TRUE, control = list(fnscale = -1))
          if(emiss_out$convergence == 0){
            subj_data[[s]]$emiss_converge[[q]][i]									<- 1
            subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ]	<- mnlHess_int(int = emiss_out$par, Obs = c(cond_y[[s]][[i]][[q]], c(1:q_emiss[q])), n_cat =  q_emiss[q])
            subj_data[[s]]$emiss_int_mle[[q]][i,]									<- emiss_out$par
          } else {
            subj_data[[s]]$emiss_converge[[q]][i]										<- 0
            subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ]	<- diag(q_emiss[q] - 1)
            subj_data[[s]]$emiss_int_mle[[q]][i,]									<- rep(0, q_emiss[q] - 1)
          }
          # if this is first iteration, use MLE for current values RW metropolis sampler
          if (iter == 2){
            emiss_c_int[[i]][[q]][s,]	<- emiss_out$par
          }
        }
      }

# Sample pouplaton values for gamma and conditional probabilities using Gibbs sampler -----------
      # gamma_mu0_n and gamma_mu_int_bar are matrices, with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
      gamma_mu0_n           <- solve(t(xx[[1]]) %*% xx[[1]] + gamma_K0)  %*% (t(xx[[1]]) %*% gamma_c_int[[i]] + gamma_K0 %*% gamma_mu0)
      gamma_V_n             <- gamma_V + t(gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) %*% (gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) + t(gamma_mu0_n - gamma_mu0) %*% gamma_K0 %*% (gamma_mu0_n - gamma_mu0)
      gamma_V_int[[i]]      <- solve(rwish(S = solve(gamma_V_n), v = gamma_nu + n_subj))
      gamma_mu_int_bar[[i]] <- gamma_mu0_n + solve(chol(t(xx[[1]]) %*% xx[[1]] + gamma_K0)) %*% matrix(rnorm((m - 1) * nx[1]), nrow = nx[1]) %*% t(solve(chol(solve(gamma_V_int[[i]]))))
      gamma_exp_int				  <- matrix(exp(c(0, gamma_mu_int_bar[[i]][1,] )), nrow  = 1)
      gamma_mu_prob_bar[[i]] 	<- gamma_exp_int / as.vector(gamma_exp_int %*% c(rep(1,(m))))

      for(q in 1:n_dep){
        emiss_mu0_n                 <- solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% (t(xx[[1 + q]]) %*% emiss_c_int[[i]][[q]] + emiss_K0[[q]] %*% emiss_mu0[[q]][[i]])
        emiss_V_n                   <- emiss_V[[q]] + t(emiss_c_int[[i]][[q]] - xx[[1 + q]] %*% emiss_mu0_n) %*% (emiss_c_int[[i]][[q]] - xx[[1 + q]] %*% emiss_mu0_n) + t(emiss_mu0_n - emiss_mu0[[q]][[i]]) %*% emiss_K0[[q]] %*% (emiss_mu0_n - emiss_mu0[[q]][[i]])
        emiss_V_int[[i]][[q]]       <- solve(rwish(S = solve(emiss_V_n), v = emiss_nu[[q]] + n_subj))
        emiss_mu_int_bar[[i]][[q]]	 <- emiss_mu0_n + solve(chol(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]])) %*% matrix(rnorm((q_emiss[q] - 1) * nx[1 + q]), nrow = nx[1 + q]) %*% t(solve(chol(solve(emiss_V_int[[i]][[q]]))))
        emiss_exp_int				       <- matrix(exp(c(0, emiss_mu_int_bar[[i]][[q]][1, ])), nrow  = 1)
        emiss_mu_prob_bar[[i]][[q]] <- emiss_exp_int / as.vector(emiss_exp_int %*% c(rep(1, (q_emiss[q]))))
      }


      # Sample subject values for gamma and conditional probabilities using RW Metropolis sampler -----------
      for (s in 1:n_subj){
        gamma_candcov_comb 			<- chol2inv(chol(subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ] + chol2inv(chol(gamma_V_int[[i]]))))
        gamma_RWout					    <- mnl_RW_once(int1 = gamma_c_int[[i]][s,], Obs = trans[[s]][[i]], n_cat = m, mu_int_bar1 = c(t(gamma_mu_int_bar[[i]]) %*% xx[[1]][s,]), V_int1 = gamma_V_int[[i]], scalar = gamma_scalar, candcov1 = gamma_candcov_comb)
        gamma[[s]][i,]  	<- PD_subj[[s]][iter, c((sum(m * q_emiss) + 1 + (i - 1) * m):(sum(m * q_emiss) + (i - 1) * m + m))] <- gamma_RWout$prob
        gamma_naccept[s, i]			<- gamma_naccept[s, i] + gamma_RWout$accept
        gamma_c_int[[i]][s,]		<- gamma_RWout$draw_int
        gamma_int_subj[[s]][iter, (1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1))] <- gamma_c_int[[i]][s,]

        start <- c(0, q_emiss * m)
        for(q in 1:n_dep){
          emiss_candcov_comb		     <- chol2inv(chol(subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ] + chol2inv(chol(emiss_V_int[[i]][[q]]))))
          emiss_RWout				       <- mnl_RW_once(int1 = emiss_c_int[[i]][[q]][s,], Obs = cond_y[[s]][[i]][[q]], n_cat = q_emiss[q], mu_int_bar1 = c(t(emiss_mu_int_bar[[i]][[q]]) %*% xx[[1 + q]][s,]), V_int1 = emiss_V_int[[i]][[q]], scalar = emiss_scalar[[q]], candcov1 = emiss_candcov_comb)
          emiss[[s]][[q]][i,]		   <- PD_subj[[s]][iter, (sum(start[1:q]) + 1 + (i - 1) * q_emiss[q]):(sum(start[1:q]) + (i - 1) * q_emiss[q] + q_emiss[q])] <- emiss_RWout$prob
          emiss_naccept[[q]][s, i]	 <- emiss_naccept[[q]][s, i] + emiss_RWout$accept
          emiss_c_int[[i]][[q]][s,] <- emiss_RWout$draw_int
          emiss_int_subj[[s]][[q]][iter, (1 + (i - 1) * (q_emiss[q] - 1)) : ((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1))]	<- emiss_c_int[[i]][[q]][s,]
        }

        if(i == m){
          delta[[s]] 		<- solve(t(diag(m) - gamma[[s]] + 1), rep(1, m))
        }
      }
    }


# End of 1 MCMC iteration, save output values --------
    gamma_int_bar[iter, ]				   	<- unlist(lapply(gamma_mu_int_bar, "[",1,))
    if(nx[1] > 1){
      gamma_cov_bar[iter, ]      	<- unlist(lapply(gamma_mu_int_bar, "[",-1,))
    }
    gamma_prob_bar[iter,]			<- unlist(gamma_mu_prob_bar)
    for(q in 1:n_dep){
      emiss_int_bar[[q]][iter, ]	<- as.vector(unlist(lapply(
        lapply(emiss_mu_int_bar, "[[", q), "[",1,)
        ))
      if(nx[1+q] > 1){
        emiss_cov_bar[[q]][iter, ]  <- as.vector(unlist(lapply(
          lapply(emiss_mu_int_bar, "[[", q), "[",-1,)
        ))
      }
      emiss_prob_bar[[q]][iter,]	<- as.vector(unlist(sapply(emiss_mu_prob_bar, "[[", q)))
    }
    if(is.whole(iter/10)){
      print(c(iter))
    }
  }


  # End of function, return output values --------
  ctime = proc.time()[3]
  print(paste("total time elapsed in minutes", round((ctime-itime) / 60, 2)))
  if(return_path == TRUE){
    out <- list(input = list(m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,
                             burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                PD_subj = PD_subj, gamma_int_subj = gamma_int_subj, emiss_int_subj = emiss_int_subj,
                gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar, emiss_int_bar = emiss_int_bar,
                emiss_cov_bar = emiss_cov_bar, gamma_prob_bar = gamma_prob_bar,
                emiss_prob_bar = emiss_prob_bar, gamma_naccept = gamma_naccept, emiss_naccept = emiss_naccept,
                sample_path = sample_path)
  } else {
    out <- list(input = list(m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,
                             burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                PD_subj = PD_subj, gamma_int_subj = gamma_int_subj, emiss_int_subj = emiss_int_subj,
                gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar, emiss_int_bar = emiss_int_bar,
                emiss_cov_bar = emiss_cov_bar, gamma_prob_bar = gamma_prob_bar,
                emiss_prob_bar = emiss_prob_bar, gamma_naccept = gamma_naccept, emiss_naccept = emiss_naccept)
  }
  class(out) <- append(class(out), "mHMM")
  return(out)
}
