#' Multilevel hidden markov model using Bayesian estimation
#'
#' This function analyses (intense longitudinal) data from multiple subjects
#' using a multilevel hidden Markov model. By using a multilevel framework, one
#' general 'population' HMM is estimated, while heterogeneity between subjects
#' is accommodated. The function can handle covariates at the subject level
#' (unlimited number), uses a hybrid metropolis within gibs sampler, and performs
#' the forward backward algorithm for all subjects in a sequential manner. Can
#' handle varying observation length over subjects

# For each parameter (also for the sub ones contained in lists), clearly specify
# the type (e.g., matrix with dimensions xx, numeric vector with lenght, string,
# list, etc)

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
#'   \item{\code{q_emis}: numeric vector with length \code{n_dep} denoting the
#'   number of observed categories for the categorical  emission distribution of
#'   each dependent variable.}}
#' @param xx List of covariates. Number of elements in the list is equal to 1 +
#'   \code{n_dep} (i.e., the number of dependent variables). The first element
#'   is used to predict the transition matrix. Subsequent elements are used to
#'   predict the emission distribution for (each of) the dependent variable(s).
#'   Each element is a matrix, with the number of rows equal to the number of
#'   subjects. The first column \emph{has to}  represent the intercept, that is,
#'   a column only consisting of ones. Subsequent columns correspond to
#'   covariates used to predict the transition matrix / emission distribution.
#'   If \code{xx} is omitted completely, \code{xx} defaults to NULL, resembling
#'   no covariates. Specific elements in the list can also be left empty (i.e.,
#'   set to NULL) to signify that either the transition probability matrix or a
#'   specific emission distribution is not predicted by covariates
#' @param start_val Start values for gamma and emis used for first run of the
#'   forward algorithm
#' @param gamma_sampler List containing start values for mle estimates of pooled
#'   data for gamma: \code{int_mle0}, \code{scalar}, and weight for the overall ll in the
#'   fractional likelihood, \code{w}.
#' @param emis_sampler List containing start values for mle estimates of pooled
#'   data for emis, emis_int_mle0, emis_scalar and weight for the overall ll in the
#'   fractional likelihood, emis_w
#' @param gamma_hyp_prior List containing \code{mu0}, \code{K0}, \code{nu} and \code{V}, see details below
#' @param emis_hyp_prior List containing \code{emis_mu0}, \code{emis_K0}, \code{emis_nu} and \code{emis_V}, see details below
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
#   \deqn{\Gamma_i ~ (idd) MNL(int_i)}
#   \deqn{int_i ~ N(mu_int_bar, V_int)}
#   \deqn{mu_int_bar ~ N(mu0, \frac{1}{K0} * V_int)}
#   \deqn{V_int ~ IW(nu, V)}
#   where \eqn{\alpha_i} denotes \code{int_i}, \eqn{\bar{\mu}_\alpha} denotes \code{mu_int_bar},  \code{} \code{} \code{} \code{} \code{}
#
#   The following hyper prior is used for the conditional emission distribution of each state at the population level:
#   \deqn{\theta_i ~ (idd) MNL(emis_int_i)}
#   \deqn{emis_int ~ N(mu_emis_int_bar, V_emis_int)}
#   \deqn{mu_emis_int_bar ~ N(emis_mu0, (1 / emis_K0) * V_emis_int)}
#   \deqn{V_emis_int 					~ IW(emis_nu, emis_V)}
#
# check thesis for correct notation and stuff
# }
#
# \subsection{Covariates}{
#   Here is some info on using covariates ... (to be completed)
# }
#
#' @return A list of ... (to be completed)
#'
#' @examples
#' # specifying general model properties
#' m <- 2
#' n_dep <- 4
#' q_emis <- c(3, 2, 3, 2)
#'
#' # specifying starting values
#' start.EM <- list(matrix(c(0.9, 0.05, 0.05, 0.05, 0.05, 0.9), byrow = TRUE,
#'                         nrow = m, ncol = q_emis[1]), # vocalizing patient
#'                  matrix(c(0.9, 0.1, 0.9, 0.1), byrow = TRUE, nrow = m,
#'                         ncol = q_emis[2]), # looking patient
#'                  matrix(c(0.05, 0.05, 0.9, 0.9, 0.05, 0.05), byrow = TRUE,
#'                         nrow = m, ncol = q_emis[3]), # vocalizing therapist
#'                  matrix(c(0.9, 0.1, 0.9, 0.1), byrow = TRUE, nrow = m,
#'                         ncol = q_emis[4])) # looking therapist
#' start.TM <- diag(.8, m)
#' start.TM[lower.tri(start.TM) | upper.tri(start.TM)] <- .2
#'
#' # run a model without covariates
#' out <- mHMM_mnl(s_data = nonverbal, gen = list(m = m, n_dep = n_dep,
#'                 q_emis = q_emis), start_val = c(as.vector(t(start.EM[[1]])),
#'                 as.vector(t(start.EM[[2]])), as.vector(t(start.EM[[3]])),
#'                 as.vector(t(start.EM[[4]])), as.vector(t(start.TM))),
#'                 mcmc = list(J = 11, burn_in = 5))
#'
#' # including covariates. Only the emission distribution for each of the 4
#' # dependent variables is predicted using standardized CDI change.
#' n_subj <- 10
#' xx <- rep(list(matrix(1, ncol = 1, nrow = n_subj)), (n_dep + 1))
#' for(i in 2:(n_dep + 1)){
#'   xx[[i]] <- cbind(xx[[i]], nonverbal_cov$std_CDI_change)
#' }
#' out2 <- mHMM_mnl(s_data = nonverbal, xx = xx, gen = list(m = m, n_dep = n_dep,
#'                 q_emis = q_emis), start_val = c(as.vector(t(start.EM[[1]])),
#'                 as.vector(t(start.EM[[2]])), as.vector(t(start.EM[[3]])),
#'                 as.vector(t(start.EM[[4]])), as.vector(t(start.TM))),
#'                 mcmc = list(J = 11, burn_in = 5))
#'
#'
# not sure if all functions given below for packages are actually still used, check!
#' @export
#' @importFrom mvtnorm dmvnorm rmvnorm dmvt rmvt
#' @importFrom MCMCpack rdirichlet rwish
#' @importFrom stats optim rnorm runif

mHMM_mnl <- function(s_data, gen, xx = NULL, start_val, gamma_sampler = NULL, emis_sampler = NULL,
                     gamma_hyp_prior = NULL, emis_hyp_prior = NULL, mcmc, return_path = FALSE){

  # Initialize data -----------------------------------
  # dependent variable(s), sample size, dimensions gamma and conditional distribuiton
  id         <- unique(s_data[,1])
  n_subj     <- length(id)
  subj_data  <- rep(list(NULL), n_subj)
  for(s in 1:n_subj){
    subj_data[[s]]$y <- s_data[s_data[,1] == id[s],][,-1]
  }
  ypooled    <- n <- NULL
  n_vary     <- numeric(n_subj)
  m          <- gen$m
  q_emis 		 <- gen$q_emis
  n_dep			 <- gen$n_dep
  emis_int_mle <- rep(list(NULL), n_dep)
  emis_mhess   <- rep(list(NULL), n_dep)
  for(q in 1:n_dep){
    emis_int_mle[[q]] <- matrix(, m, (q_emis[q] - 1))
    emis_mhess[[q]] <- matrix(, (q_emis[q] - 1) * m, (q_emis[q] - 1))
  }
  for(s in 1:n_subj){
    ypooled   <- rbind(ypooled, subj_data[[s]]$y)
    n         <- dim(subj_data[[s]]$y)[1]
    n_vary[s] <- n
    subj_data[[s]]	<- c(subj_data[[s]], n = n, list(converge = numeric(m), int_mle = matrix(, m, (m - 1)),
                                                    mhess = matrix(, (m - 1) * m, (m - 1)), emis_converge =
                                                      rep(list(numeric(m)), n_dep), emis_int_mle = emis_int_mle, emis_mhess = emis_mhess))
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
      }
    }
  }

  # Initialize mcmc argumetns
  J 				<- mcmc$J
  burn_in			<- mcmc$burn_in


  # Initalize priors and hyper priors --------------------------------
  # Initialize gamma sampler
  if(is.null(gamma_sampler)) {
    int_mle0  <- mu_prop <- rep(0, m - 1)
    scalar    <- 2.93 / sqrt(m - 1)
    w         <- .1
  } else {
    int_mle0  <- gamma_sampler$int_mle0
    mu_prop   <- gamma_sampler$mu_prop
    scalar    <- gamma_sampler$scalar
    w         <- gamma_sampler$w
  }

  # Initialize emis sampler
  if(is.null(emis_sampler)){
    emis_int_mle0 <- rep(list(NULL), n_dep)
    emis_scalar	<- rep(list(NULL), n_dep)
    mu_emis_prop 	<-  rep(list(NULL), n_dep)
    for(q in 1:n_dep){
      emis_int_mle0[[q]]	<- mu_emis_prop[[q]] <- rep(0, q_emis[q] - 1)
      emis_scalar[[q]] 	<- 2.93 / sqrt(q_emis[q] - 1)
    }
    emis_w		<- .1
  } else {
    emis_int_mle0	<- emis_sampler$emis_int_mle0
    mu_emis_prop 	<- emis_sampler$mu_emis_prop
    emis_scalar 	<- emis_sampler$emis_scalar
    emis_w    		<- emis_sampler$emis_w
  }

  # Initialize Gamma hyper prior
  if(is.null(gamma_hyp_prior)){
    mu0			<- matrix(0,nrow = nx, ncol = m - 1)
    K0			<- diag(1, nx[1])
    nu			<- 3 + m - 1
    V			  <- nu * diag(m - 1)
  } else {
    ###### BUILD in a warning / check if mu0 is a matrix when given, with  nrows equal to the number of covariates
    mu0			<- gamma_hyp_prior$mu0
    K0			<- gamma_hyp_prior$K0
    nu			<- gamma_hyp_prior$nu
    V			  <- gamma_hyp_prior$V
  }

  # Initialize Pr hyper prior
  # emis_mu0: for each dependent variable, emis_mu0 is a list, with one element for each state.
  # Each element is a matrix, with number of rows equal to the number of covariates (with the intercept being one cov),
  # and the number of columns equal to q_emis[q] - 1.
  emis_mu0	  <- rep(list(vector("list", m)), n_dep)
  emis_nu	    <- rep(list(NULL), n_dep)
  emis_V	    <- rep(list(NULL), n_dep)
  emis_K0     <- rep(list(NULL), n_dep)
  if(is.null(emis_hyp_prior)){
    for(q in 1:n_dep){
      for(i in 1:m){
        emis_mu0[[q]][[i]]		<- matrix(0, ncol = q_emis[q] - 1, nrow = nx[1 + q])
      }
      emis_nu[[q]]		<- 3 + q_emis[q] - 1
      emis_V[[q]]		  <- emis_nu[[q]] * diag(q_emis[q] - 1)
      emis_K0[[q]]		<- diag(1, nx[1 + q])
    }
  } else {
    for(q in 1:n_dep){
      # emis_hyp_prior[[q]]$emis_mu0 has to contain a list with lenght equal to m, and each list contains matrix with number of rows equal to number of covariates for that dep. var.
      # stil build in a CHECK, with warning / stop / switch to default prior
      emis_mu0[[q]]	 <- emis_hyp_prior[[q]]$emis_mu0
      emis_nu[[q]]	 <- emis_hyp_prior[[q]]$emis_nu
      emis_V[[q]]		 <- emis_hyp_prior[[q]]$emis_V
      emis_K0[[q]]	 <- diag(emis_hyp_prior$emis_K0, nx[1 + q])
    }
  }


  # Define objects used to store data in mcmc algorithm, not returned ----------------------------
  # overall
  c <- llk <- numeric(1)
  sample_path <- lapply(n_vary, dif_matrix, cols = J)
  trans <- rep(list(vector("list", m)), n_subj)

  # gamma
  mle_pooled <-int_mle_pooled <- mhess_pooled <- pooled_ll <- vector("list", m)
  c_int <- rep(list(matrix(, n_subj, (m-1))), m)
  mu_int_bar <- V_int <- vector("list", m)
  mu_gamma_prob_bar <- rep(list(numeric(m)), m)
  naccept <- matrix(0, n_subj, m)

  # emis
  cond_y <- lapply(rep(n_dep, n_subj), nested_list, m = m)
  emis_mle_pooled <- emis_int_mle_pooled <- emis_mhess_pooled <- emis_pooled_ll <- rep(list(vector("list", n_dep)), m)
  emis_c_int <- rep(list(lapply(q_emis - 1, dif_matrix, rows = n_subj)), m)
  mu_emis_int_bar <- V_emis_int <- rep(list(vector("list", n_dep)), m)
  mu_emis_prob_bar <- rep(list(lapply(q_emis, dif_vector)), m)
  emis_naccept <- rep(list(matrix(0, n_subj, m)), n_dep)


  # Define objects that are returned from mcmc algorithm ----------------------------
  # Define object for subject specific posterior density, put start values on first row
  PD 					  <- matrix(, nrow = J, ncol = sum(m * q_emis) + m * m + 1)
  PD_emis_names   <- paste("q", 1, "_emis", rep(1:q_emis[1], m), "_S", rep(1:m, each = q_emis[1]), sep = "")
  for(q in 2:n_dep){
    PD_emis_names <- c(PD_emis_names, paste("q", q, "_emis", rep(1:q_emis[q], m), "_S", rep(1:m, each = q_emis[q]), sep = ""))
  }
  colnames(PD) 	<- c(PD_emis_names, paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = ""), "LL")
  PD[1, 1:((sum(m * q_emis) + m * m))] <- start_val[1:((sum(m * q_emis) + m * m))]
  PD_subj				<- rep(list(PD), n_subj)

  # Define object for subject specific posterior density (regression coefficients parameterization )
  int_subj			<- rep(list(matrix(, nrow = J, ncol = (m-1) * m)), n_subj)
  emis_int_subj			<- rep(list(lapply((q_emis - 1) * m, dif_matrix, rows = J)), n_subj)

  # Define object for population posterior density ('normal' and regression coefficients parameterization )
  gamma_prob_bar		<- matrix(, nrow = J, ncol = (m * m))
  emis_prob_bar			<- lapply(q_emis * m, dif_matrix, rows = J)
  int_bar				<- matrix(, nrow = J, ncol = ((m-1) * m) * nx[1])
  emis_int_bar			<- lapply((q_emis-1) * m * nx[-1], dif_matrix, rows = J)

  # Put starting values in place for fist run forward algorithm
  emis_sep 			<- vector("list", n_dep)
  for(q in 1:n_dep){
    start <- c(0, q_emis * m)
    emis_sep[[q]] <- matrix(start_val[(sum(start[1:q]) + 1):(sum(start[1:q]) + (m * q_emis[q]))], byrow = TRUE, ncol = q_emis[q], nrow = m)
  }
  emis				  <- rep(list(emis_sep), n_subj)
  gamma 			<- rep(list(matrix(start_val[(sum(m*q_emis) + 1):(sum(m * q_emis) + m * m)], byrow = TRUE, ncol = m)), n_subj)
  delta 			<- rep(list(solve(t(diag(m) - gamma[[1]] + 1), rep(1, m))), n_subj)


  # Start analysis --------------------------------------------
  # Run the MCMC algorithm
  itime <- proc.time()[3]
  for (iter in 2 : J){

    # For each subject, obtain sampled state sequence with subject individual parameters ----------
    for(s in 1:n_subj){
      # Run forward algorithm, obtain subject specific forward proababilities and log likelihood
      forward				<- cat_Mult_HMM_fw(x = subj_data[[s]]$y, m = m, emis = emis[[s]], gamma = gamma[[s]], n_dep = n_dep, delta=NULL)
      alpha         <- forward$forward_p
      c             <- max(forward$la[, subj_data[[s]]$n])
      llk           <- c + log(sum(exp(forward$la[, subj_data[[s]]$n] - c)))
      PD_subj[[s]][iter, sum(m * q_emis) + m * m + 1] <- llk

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
          cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q], 1:q_emis[q])
        }
      }
    }

    # The remainder of the mcmc algorithm is state specific
    for(i in 1:m){

      # Obtain MLE of the covariance matrices and log likelihood of gamma and emis at subject and population level -----------------
      # used to scale the propasal distribution of the RW Metropolis sampler

      # population level, transition matrix
      trans_pooled			  <- factor(c(unlist(sapply(trans, "[[", i)), c(1:m)))
      mle_pooled[[i]] 		<- optim(int_mle0, llmnl_int, Obs = trans_pooled, n_cat = m, method = "BFGS", hessian = TRUE, control = list(fnscale = -1))
      int_mle_pooled[[i]] <- mle_pooled[[i]]$par
      pooled_ll[[i]]			<- mle_pooled[[i]]$value

      # population level, conditional probabilities, seperate for each dependent variable
      for(q in 1:n_dep){
        cond_y_pooled					      <- numeric()
        ### MOET OOK ECHT BETER KUNNEN, eerst # cond_y_pooled				<- unlist(sapply(cond_y, "[[", m))
        for(s in 1:n_subj){
          cond_y_pooled             <- c(cond_y_pooled, cond_y[[s]][[i]][[q]])
        }
        emis_mle_pooled[[i]][[q]]		  <- optim(emis_int_mle0[[q]], llmnl_int, Obs = c(cond_y_pooled, c(1:q_emis[q])), n_cat = q_emis[q],
                                              method = "BFGS", hessian = TRUE, control = list(fnscale = -1))
        emis_int_mle_pooled[[i]][[q]]	<- emis_mle_pooled[[i]][[q]]$par
        emis_pooled_ll[[i]][[q]]			<- emis_mle_pooled[[i]][[q]]$value
      }

      # subject level, transition matrix
      for (s in 1:n_subj){
        wgt 				<- subj_data[[s]]$n / n_total
        out					<- optim(int_mle_pooled[[i]], llmnl_int_frac, Obs = c(trans[[s]][[i]], c(1:m)), n_cat = m,
                         pooled_likel = pooled_ll[[i]], w = w, wgt = wgt, method="BFGS", hessian = TRUE, control = list(fnscale = -1))
        if(out$convergence == 0){
          subj_data[[s]]$converge[i]									<- 1
          subj_data[[s]]$mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ]	<- mnlHess_int(int = out$par, Obs = c(trans[[s]][[i]], c(1:m)), n_cat =  m)
          subj_data[[s]]$int_mle[i,]									<- out$par
        } else {
          subj_data[[s]]$converge[i]									<- 0
          subj_data[[s]]$mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)),]	<- diag(m-1)
          subj_data[[s]]$int_mle[i,]									<- rep(0, m - 1)
        }
        # if this is first iteration, use MLE for current values RW metropolis sampler
        if (iter == 2){
          c_int[[i]][s,]		<- out$par
        }

        # subject level, conditional probabilities, seperate for each dependent variable
        for(q in 1:n_dep){
          emis_out				<- optim(emis_int_mle_pooled[[i]][[q]], llmnl_int_frac, Obs = c(cond_y[[s]][[i]][[q]], c(1:q_emis[q])),
                               n_cat = q_emis[q], pooled_likel = emis_pooled_ll[[i]][[q]], w = emis_w, wgt = wgt, method = "BFGS", hessian = TRUE, control = list(fnscale = -1))
          if(emis_out$convergence == 0){
            subj_data[[s]]$emis_converge[[q]][i]									<- 1
            subj_data[[s]]$emis_mhess[[q]][(1 + (i - 1) * (q_emis[q] - 1)):((q_emis[q] - 1) + (i - 1) * (q_emis[q] - 1)), ]	<- mnlHess_int(int = emis_out$par, Obs = c(cond_y[[s]][[i]][[q]], c(1:q_emis[q])), n_cat =  q_emis[q])
            subj_data[[s]]$emis_int_mle[[q]][i,]									<- emis_out$par
          } else {
            subj_data[[s]]$emis_converge[[q]][i]										<- 0
            subj_data[[s]]$emis_mhess[[q]][(1 + (i - 1) * (q_emis[q] - 1)):((q_emis[q] - 1) + (i - 1) * (q_emis[q] - 1)), ]	<- diag(q_emis[q] - 1)
            subj_data[[s]]$emis_int_mle[[q]][i,]									<- rep(0, q_emis[q] - 1)
          }
          # if this is first iteration, use MLE for current values RW metropolis sampler
          if (iter == 2){
            emis_c_int[[i]][[q]][s,]	<- emis_out$par
          }
        }
      }

      # Sample pouplaton values for gamma and conditional probabilities using Gibbs sampler -----------
      # mu0_n and mu_int_bar are matrices, with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
      mu0_n           <- solve(t(xx[[1]]) %*% xx[[1]] + K0)  %*% (t(xx[[1]]) %*% c_int[[i]] + K0 %*% mu0)
      V_n             <- V + t(c_int[[i]] - xx[[1]] %*% mu0_n) %*% (c_int[[i]] - xx[[1]] %*% mu0_n) + t(mu0_n - mu0) %*% K0 %*% (mu0_n - mu0)
      V_int[[i]]      <- solve(rwish(S = solve(V_n), v = nu + n_subj))
      mu_int_bar[[i]] <- mu0_n + solve(chol(t(xx[[1]]) %*% xx[[1]] + K0)) %*% matrix(rnorm((m - 1) * nx[1]), nrow = nx[1]) %*% t(solve(chol(solve(V_int[[i]]))))
      exp_int				  <- matrix(exp(c(0, mu_int_bar[[i]][1,] )), nrow  = 1)
      mu_gamma_prob_bar[[i]] 	<- exp_int / as.vector(exp_int %*% c(rep(1,(m))))

      for(q in 1:n_dep){
        emis_mu0_n                 <- solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emis_K0[[q]]) %*% (t(xx[[1 + q]]) %*% emis_c_int[[i]][[q]] + emis_K0[[q]] %*% emis_mu0[[q]][[i]])
        emis_V_n                   <- emis_V[[q]] + t(emis_c_int[[i]][[q]] - xx[[1 + q]] %*% emis_mu0_n) %*% (emis_c_int[[i]][[q]] - xx[[1 + q]] %*% emis_mu0_n) + t(emis_mu0_n - emis_mu0[[q]][[i]]) %*% emis_K0[[q]] %*% (emis_mu0_n - emis_mu0[[q]][[i]])
        V_emis_int[[i]][[q]]       <- solve(rwish(S = solve(emis_V_n), v = emis_nu[[q]] + n_subj))
        mu_emis_int_bar[[i]][[q]]	 <- emis_mu0_n + solve(chol(t(xx[[1 + q]]) %*% xx[[1 + q]] + emis_K0[[q]])) %*% matrix(rnorm((q_emis[q] - 1) * nx[1 + q]), nrow = nx[1 + q]) %*% t(solve(chol(solve(V_emis_int[[i]][[q]]))))
        exp_emis_int				       <- matrix(exp(c(0, mu_emis_int_bar[[i]][[q]][1, ])), nrow  = 1)
        mu_emis_prob_bar[[i]][[q]] <- exp_emis_int / as.vector(exp_emis_int %*% c(rep(1, (q_emis[q]))))
      }


      # Sample subject values for gamma and conditional probabilities using RW Metropolis sampler -----------
      for (s in 1:n_subj){
        candcov_comb 			<- chol2inv(chol(subj_data[[s]]$mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ] + chol2inv(chol(V_int[[i]]))))
        RWout					    <- mnl_RW_once(int1 = c_int[[i]][s,], Obs = trans[[s]][[i]], n_cat = m, mu_int_bar1 = c(t(mu_int_bar[[i]]) %*% xx[[1]][s,]), V_int1 = V_int[[i]], scalar = scalar, candcov1 = candcov_comb)
        gamma[[s]][i,]  	<- PD_subj[[s]][iter, c((sum(m * q_emis) + 1 + (i - 1) * m):(sum(m * q_emis) + (i - 1) * m + m))] <- RWout$prob
        naccept[s, i]			<- naccept[s, i] + RWout$accept
        c_int[[i]][s,]		<- RWout$draw_int
        int_subj[[s]][iter, (1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1))] <- c_int[[i]][s,]

        start <- c(0, q_emis * m)
        for(q in 1:n_dep){
          emis_candcov_comb		     <- chol2inv(chol(subj_data[[s]]$emis_mhess[[q]][(1 + (i - 1) * (q_emis[q] - 1)):((q_emis[q] - 1) + (i - 1) * (q_emis[q] - 1)), ] + chol2inv(chol(V_emis_int[[i]][[q]]))))
          emis_RWout				       <- mnl_RW_once(int1 = emis_c_int[[i]][[q]][s,], Obs = cond_y[[s]][[i]][[q]], n_cat = q_emis[q], mu_int_bar1 = c(t(mu_emis_int_bar[[i]][[q]]) %*% xx[[1 + q]][s,]), V_int1 = V_emis_int[[i]][[q]], scalar = emis_scalar[[q]], candcov1 = emis_candcov_comb)
          emis[[s]][[q]][i,]		   <- PD_subj[[s]][iter, (sum(start[1:q]) + 1 + (i - 1) * q_emis[q]):(sum(start[1:q]) + (i - 1) * q_emis[q] + q_emis[q])] <- emis_RWout$prob
          emis_naccept[[q]][s, i]	 <- emis_naccept[[q]][s, i] + emis_RWout$accept
          emis_c_int[[i]][[q]][s,] <- emis_RWout$draw_int
          emis_int_subj[[s]][[q]][iter, (1 + (i - 1) * (q_emis[q] - 1)) : ((q_emis[q] - 1) + (i - 1) * (q_emis[q] - 1))]	<- emis_c_int[[i]][[q]][s,]
        }

        if(i == m){
          delta[[s]] 		<- solve(t(diag(m) - gamma[[s]] + 1), rep(1, m))
        }
      }
    }


    # End of 1 MCMC iteration, save output values --------
    int_bar[iter, ]				   	<- unlist(mu_int_bar)
    gamma_prob_bar[iter,]			<- unlist(mu_gamma_prob_bar)
    for(q in 1:n_dep){
      emis_int_bar[[q]][iter, ]	<- as.vector(unlist(sapply(mu_emis_int_bar, "[[", q)))
      emis_prob_bar[[q]][iter,]	<- as.vector(unlist(sapply(mu_emis_prob_bar, "[[", q)))
    }
    if(is.whole(iter/10)){
      print(c(iter))
    }
  }


  # End of function, return output values --------
  ctime = proc.time()[3]
  print(paste("total time elapsed in minutes", round((ctime-itime) / 60, 2)))
  if(return_path == TRUE){
    return(list(PD_subj = PD_subj, int_subj = int_subj, emis_int_subj = emis_int_subj, int_bar = int_bar, emis_int_bar = emis_int_bar, gamma_prob_bar = gamma_prob_bar, emis_prob_bar = emis_prob_bar, naccept = naccept, emis_naccept = emis_naccept, sample_path = sample_path))
  } else {
    return(list(PD_subj = PD_subj, int_subj = int_subj, emis_int_subj = emis_int_subj, int_bar = int_bar, emis_int_bar = emis_int_bar, gamma_prob_bar = gamma_prob_bar, emis_prob_bar = emis_prob_bar, naccept = naccept, emis_naccept = emis_naccept))
  }
}
