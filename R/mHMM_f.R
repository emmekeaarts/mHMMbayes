#' Multilevel hidden  Markov model using Bayesian estimation
#'
#' \code{mHMM_f} fits a multilevel (also known as mixed or random effects)
#' hidden Markov model (HMM) to intense longitudinal data with (ONLY) continuous
#' (i.e., normally distributed), observations of multiple subjects using
#' Bayesian estimation, and creates an object of class \code{mHMM_f}. By using a
#' multilevel framework, we allow for heterogeneity in the transition model
#' parameters between subjects, while estimating one overall HMM. The function
#' includes the possibility to add covariates at level 2 (i.e., at the subject
#' level) and have varying observation lengths over subjects. For a short
#' description of the package see \link{mHMMbayes}. See
#' \code{vignette("tutorial-mhmm")} for an introduction to multilevel hidden
#' Markov models and the package, and see \code{vignette("estimation-mhmm")} for
#' an overview of the used estimation algorithms.
#'
#' Covariates specified in \code{xx} can either be dichotomous or continuous
#' variables. Dichotomous variables have to be coded as 0/1 variables.
#' Categorical or factor variables can as yet not be used as predictor
#' covariates. The user can however break up the categorical variable in
#' multiple dummy variables (i.e., dichotomous variables), which can be used
#' simultaneously in the analysis. Continuous predictors are automatically
#' centered. That is, the mean value of the covariate is subtracted from all
#' values of the covariate such that the new mean equals zero. This is done such
#' that the presented probabilities in the output (i.e., for the population
#' transition probability matrix and population emission probabilities)
#' correspond to the predicted probabilities at the average value of the
#' covariate(s).
#'
#'
#' @param s_data A matrix containing the observations to be modeled, where the
#'   rows represent the observations over time. In \code{s_data}, the first
#'   column indicates subject id number. Hence, the id number is repeated over
#'   rows equal to the number of observations for that subject. The subsequent
#'   columns contain the dependent variable(s). The total number of rows are
#'   equal to the sum over the number of observations of each subject, and the
#'   number of columns are equal to the number of dependent variables
#'   (\code{n_dep}) + 1. The number of observations can vary over subjects.
#' @param data_distr String vector with length 1 describing the
#'   observation type of the data. Currently supported are \code{'categorical'}
#'   , \code{'continuous'}, and \code{'count'}. Note that for multivariate
#'   data, all dependent variables are assumed to be of the same observation
#'   type. The default equals to \code{data_distr = 'categorical'}.
#' @param gen List containing the following elements denoting the general model
#'   properties:
#'   \itemize{
#'   \item{\code{m}: numeric vector with length 1 denoting the number
#'   of hidden states}
#'   \item{\code{n_dep}: numeric vector with length 1 denoting the
#'   number of dependent variables}}
#' @param xx An optional list of (level 2) covariates to predict the transition
#'   matrix. Level 2 covariate(s) means that there is one observation per
#'   subject of each covariate. The first element in the list \code{xx} is used
#'   to predict the transition matrix. Each element in the list is a matrix,
#'   with the number of rows equal to the number of subjects. The first column
#'   of each matrix represents the intercept, that is, a column only consisting
#'   of ones. Subsequent columns correspond to covariates used to predict the
#'   transition matrix / emission distribution. See \emph{Details} for more
#'   information on the use of covariates. Note that since the emission means
#'   are fixed over subjects in this version, no covariates can be used to
#'   predict variation in the emission means over subjects. Hence, this list
#'   here always contains only 1 element.
#'
#'   If \code{xx} is omitted completely, \code{xx} defaults to \code{NULL},
#'   resembling no covariates.
#' @param start_val List containing the start values for the transition
#'   probability matrix gamma and the emission distribution(s). The first
#'   element of the list contains a \code{m} by \code{m} matrix with the start
#'   values for gamma. The subsequent \code{n_dep} elements each contain a
#'   matrix with the start values for the emission distribution(s): a \code{m}
#'   by 2 matrix denoting the mean (first column) and standard deviation (second
#'   column) of the Normal emission distribution within each state (rows). Note
#'   that \code{start_val} should not contain nested lists (i.e., lists within
#'   lists).
#' @param mcmc List of Markov chain Monte Carlo (MCMC) arguments, containing the
#'   following elements:
#'   \itemize{\item{\code{J}: numeric vector with length 1 denoting the number
#'   of iterations of the MCMC algorithm}
#'   \item{\code{burn_in}: numeric vector with length 1 denoting the
#'   burn-in period for the MCMC algorithm.}}
#' @param return_path A logical scalar. Should the sampled state sequence
#'   obtained at each iteration and for each subject be returned by the function
#'   (\code{sample_path = TRUE}) or not (\code{sample_path = FALSE}). Note that
#'   the sampled state sequence is quite a large object, hence the default
#'   setting is \code{sample_path = FALSE}. Can be used for local decoding
#'   purposes.
#' @param show_progress A logical scaler. Should the function show a text
#'   progress bar in the \code{R} console to represent the progress of the
#'   algorithm (\code{show_progress = TRUE}) or not (\code{show_progress =
#'   FALSE}). Defaults to \code{show_progress = TRUE}.
#' @param gamma_hyp_prior An optional object of class \code{mHMM_prior_gamma}
#'   containing user specified parameter values for the hyper-prior distribution
#'   on the transition probability matrix gamma, generated by the function
#'   \code{\link{prior_gamma}}.
#' @param emiss_hyp_prior An object of the class \code{mHMM_prior_emiss}
#'   containing user specified parameter values for the hyper-prior distribution
#'   on the (fixed over subjects) Normal emission distribution(s), generated by
#'   the function \code{\link{prior_emiss_cont_f}}.
#' @param gamma_sampler An optional object of the class \code{mHMM_pdRW_gamma}
#'   containing user specified settings for the proposal distribution of the
#'   random walk (RW) Metropolis sampler on the subject level transition
#'   probability matrix parameters, generated by the function
#'   \code{\link{pd_RW_gamma}}.
#'
#' @return \code{mHMM_f} returns an object of class \code{mHMM_f}, which has
#'   \code{print} and \code{summary} methods to see the results.
#'   The object contains the following components:
#'   \describe{
#'   \item{\code{PD_subj}}{A list containing one list per subject with the
#'   elements \code{trans_prob} and \code{log_likl}, providing the subject
#'   parameter estimates over the iterations of the MCMC sampler.
#'   \code{trans_prob} relates to the transition probabilities gamma,
#'   and \code{log_likl} to the log likelihood over the MCMC iterations.
#'   Iterations are contained in the rows, the parameters in the columns.}
#'   \item{\code{gamma_prob_bar}}{A matrix containing the group level parameter
#'   estimates of the transition probabilities over the iterations of the hybrid
#'   Metropolis within Gibbs sampler. The iterations of the sampler are
#'   contained in the rows, and the columns contain the group level parameter
#'   estimates. If covariates were included in the analysis, the group level
#'   probabilities represent the predicted probability given that the covariate
#'   is at the average value for continuous covariates, or given that the
#'   covariate equals zero for dichotomous covariates.}
#'   \item{\code{gamma_int_bar}}{A matrix containing the group level intercepts
#'   of the Multinomial logistic regression modeling the transition
#'   probabilities over the iterations of the hybrid Metropolis within Gibbs
#'   sampler. The iterations of the sampler are contained in the rows, and the
#'   columns contain the group level intercepts.}
#'   \item{\code{gamma_cov_bar}}{A matrix containing the group level regression
#'   coefficients of the Multinomial logistic regression predicting the
#'   transition probabilities over the iterations of the hybrid Metropolis within
#'   Gibbs sampler. The iterations of the sampler are contained in the rows, and
#'   the columns contain the group level regression coefficients.}
#'   \item{\code{gamma_V_int_bar}}{A matrix containing the variance
#'   components for the subject-level intercepts (between subject variances)
#'   of the multinomial logistic regression modeling the transition
#'   probabilities over the iterations of the hybrid Metropolis within Gibbs
#'   sampler. The iterations of the sampler are contained in the rows, and the
#'   columns contain the variance components for the subject level intercepts.
#'   Note that only the intercept variances (and not the co-variances) are
#'   returned.}
#'   \item{\code{gamma_int_subj}}{A list containing one matrix per subject
#'   denoting the subject level intercepts of the Multinomial logistic
#'   regression modeling the transition probabilities over the iterations of the
#'   hybrid Metropolis within Gibbs sampler. The iterations of the sampler are
#'   contained in the rows, and the columns contain the subject level
#'   intercepts.}
#'   \item{\code{gamma_naccept}}{A matrix containing the number of accepted
#'   draws at the subject level RW Metropolis step for each set of parameters of
#'   the transition probabilities. The subjects are contained in the rows, and
#'   the columns contain the sets of parameters.}
#'   \item{\code{emiss_mu}}{A list containing one matrix per dependent
#'   variable, denoting the (fixed over subjects) means of the Normal emission
#'   distribution of each dependent variable over the iterations of the Gibbs
#'   sampler. The iterations of the sampler are contained in the rows of the
#'   matrix, and the columns contain the group level emission means.}
#'   \item{\code{emiss_sd}}{A list containing one matrix per dependent
#'   variable, denoting the (fixed over subjects) standard deviation of the
#'   Normal emission distributions over the iterations of the Gibbs sampler. The
#'   iterations of the sampler are contained in the rows of the matrix, and the
#'   columns contain the group level emission variances.}
#'   \item{\code{input}}{Overview of used input specifications: the number of
#'   states \code{m}, the number of used dependent variables \code{n_dep}, the
#'   number of iterations \code{J} and the specified burn in period
#'   \code{burn_in} of the hybrid Metropolis within Gibbs sampler, the number of
#'   subjects \code{n_subj}, the observation length for each subject
#'   \code{n_vary}, and the column names of the dependent variables
#'   \code{dep_labels}.}
#'   \item{\code{input}}{Overview of used input specifications: the distribution
#'   type of the observations \code{data_distr}, the number of
#'   states \code{m}, the number of used dependent variables \code{n_dep}, the number of iterations
#'   \code{J} and the specified burn in period \code{burn_in} of the hybrid
#'   Metropolis within Gibbs sampler, the number of subjects \code{n_subj}, the
#'   observation length for each subject \code{n_vary}, and the column names of
#'   the dependent variables \code{dep_labels}.}
#'   \item{\code{sample_path}}{A list containing one matrix per subject with the
#'   sampled hidden state sequence over the hybrid Metropolis within Gibbs
#'   sampler. The time points of the dataset are contained in the rows, and the
#'   sampled paths over the iterations are contained in the columns. Only
#'   returned if \code{return_path = TRUE}. }
#'   }
#'
#'
#' @seealso \code{\link{sim_mHMM}} for simulating multilevel hidden Markov data,
#'   \code{\link{vit_mHMM}} for obtaining the most likely hidden state sequence
#'   for each subject using the Viterbi algorithm, \code{\link{obtain_gamma}}
#'   and \code{\link{obtain_emiss}} for obtaining the transition or emission
#'   distribution probabilities of a fitted model at the group or subject level,
#'   and \code{\link{plot.mHMM}} for plotting the posterior densities of a
#'   fitted model.
#'
#' @references
#' \insertRef{rabiner1989}{mHMMbayes}
#'
#' \insertRef{scott2002}{mHMMbayes}
#'
#' \insertRef{altman2007}{mHMMbayes}
#'
#' \insertRef{rossi2012}{mHMMbayes}
#'
#' \insertRef{zucchini2017}{mHMMbayes}
#'
#' @examples
#' ###### Example on continuous simulated data
#' \donttest{ # simulating multivariate continuous data
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
#' data_cont <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous',
#'                       gen = list(m = m, n_dep = n_dep),
#'                       gamma = gamma, emiss_distr = emiss_distr,
#'                       var_gamma = .1, var_emiss = c(5^2, 0.2^2))
#'
#' # Specify hyper-prior for the continuous emission distribution
#' manual_prior_emiss <- prior_emiss_cont_f(
#'                         gen = list(m = m, n_dep = n_dep),
#'                         emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
#'                                          matrix(c(7, 8, 18), nrow = 1)),
#'                         emiss_K0 = list(1, 1),
#'                         emiss_V =  list(rep(5^2, m), rep(0.5^2, m)),
#'                         emiss_nu = list(1, 1))
#'
#' # Run the model on the simulated data:
#' # Note that for reasons of running time, J is set at a ridiculous low value.
#' # One would typically use a number of iterations J of at least 1000,
#' # and a burn_in of 200.
#' out_3st_cont_sim <- mHMM_f(s_data = data_cont$obs,
#'                          data_distr = 'continuous',
#'                          gen = list(m = m, n_dep = n_dep),
#'                          start_val = c(list(gamma), emiss_distr),
#'                          emiss_hyp_prior = manual_prior_emiss,
#'                          mcmc = list(J = 11, burn_in = 5))
#'
#' }
#'
#'
#' @export
#'
#'

mHMM_f <- function(s_data, data_distr = 'continuous', gen, xx = NULL, start_val, mcmc, return_path = FALSE, show_progress = TRUE,
                 gamma_hyp_prior = NULL, emiss_hyp_prior = NULL, gamma_sampler = NULL){
  # Initialize data -----------------------------------
  # dependent variable(s), sample size, dimensions gamma and conditional distribution
  if(sum(objects(gen) %in% "m") != 1 | sum(objects(gen) %in% "n_dep") != 1){
    stop("The input argument gen should contain the elements m and n_dep.")
  }
  if(data_distr != 'continuous'){
    stop("this developers version can only be used for observations that have a Normal (i.e., Gaussian) emission distribution")
  }
  n_dep			 <- gen$n_dep
  dep_labels <- colnames(s_data[,2:(n_dep+1)])
  id         <- unique(s_data[,1])
  n_subj     <- length(id)
  subj_data  <- rep(list(NULL), n_subj)
  if(sum(sapply(s_data, is.factor)) > 0 ){
    stop("Your data contains factorial variables, which cannot be used as input in the function mHMM. All variables have to be numerical.")
  }
  if(sum(class(s_data) %in% "list") > 0){
    stop("The input data specified in s_data should be a matrix or data frame")
  }
  for(s in 1:n_subj){
    subj_data[[s]]$y <- as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep)
  }
  ypooled    <- n_t <- NULL
  n_vary     <- numeric(n_subj)
  m          <- gen$m
  for(s in 1:n_subj){
    ypooled   <- rbind(ypooled, subj_data[[s]]$y)
    n_t       <- dim(subj_data[[s]]$y)[1]
    n_vary[s] <- n_t
    subj_data[[s]]	<- c(subj_data[[s]], n_t = n_t, list(gamma_mhess = matrix(NA_real_, (m - 1) * m, (m - 1))))
  }
  n_total 		<- dim(ypooled)[1]

  # covariates
  n_dep1 <- 1 + n_dep
  nx <- numeric(n_dep1)
  if (is.null(xx)){
    xx <- rep(list(matrix(1, ncol = 1, nrow = n_subj)), n_dep1)
    nx[] <- 1
  } else {
    if(!is.list(xx) | length(xx) != n_dep1){
      stop("If xx is specified, xx should be a list, with the number of elements equal to the number of dependent variables + 1")
    }
    for(i in 1:n_dep1){
      if (is.null(xx[[i]])){
        xx[[i]] <- matrix(1, ncol = 1, nrow = n_subj)
        nx[i] <- 1
      } else {
        nx[i] <- ncol(xx[[i]])
        if(sum(is.na(xx[[i]])) > 1){
          stop("Currently, missing values (NA) are only allowed in the dependent variables and not in the covariate(s) values xx")
        }
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
  burn_in		<- mcmc$burn_in


  # Initalize priors and hyper priors --------------------------------
  # Initialize gamma sampler
  if(m > 1){
    if(is.null(gamma_sampler)) {
      gamma_int_mle0  <- matrix(0, nrow = m, ncol = m - 1)
      gamma_scalar    <- 2.93 / sqrt(m - 1)
      gamma_w         <- .1
    } else {
      if (!is.mHMM_pdRW_gamma(gamma_sampler)){
        stop("The input object specified for gamma_sampler should be from the class mHMM_pdRW_gamma, obtained by using the function pd_RW_gamma")
      }
      if (gamma_sampler$m != m){
        stop("The number of states specified in m is not equal to the number of states specified when setting the proposal distribution of the RW Metropolis sampler on gamma using the function pd_RW_gamma")
      }
      gamma_int_mle0  <- gamma_sampler$gamma_int_mle0
      gamma_scalar    <- gamma_sampler$gamma_scalar
      gamma_w         <- gamma_sampler$gamma_w
    }
  }

  # Initialize Gamma hyper prior
  if(m > 1){
    if(is.null(gamma_hyp_prior)){
      gamma_mu0	  <- rep(list(matrix(0,nrow = nx[1], ncol = m - 1)), m)
      gamma_K0			<- diag(1, nx[1])
      gamma_nu			<- 3 + m - 1
      gamma_V			  <- gamma_nu * diag(m - 1)
    } else {
      if (!is.mHMM_prior_gamma(gamma_hyp_prior)){
        stop("The input object specified for gamma_hyp_prior should be from the class mHMM_prior_gamma, obtained by using the function prior_gamma.")
      }
      if (gamma_hyp_prior$m != m){
        stop("The number of states specified in m is not equal to the number of states specified when creating the informative hper-prior distribution gamma using the function prior_gamma.")
      }
      if(is.null(gamma_hyp_prior$n_xx_gamma) & nx[1] > 1){
        stop("Covariates were specified to predict gamma, but no covariates were specified when creating the informative hyper-prior distribution on gamma using the function prior_gamma.")
      }
      if(!is.null(gamma_hyp_prior$n_xx_gamma)){
        if(gamma_hyp_prior$n_xx_gamma != (nx[1] - 1)){
          stop("The number of covariates specified to predict gamma is not equal to the number of covariates specified when creating the informative hper-prior distribution on gamma using the function prior_gamma.")
        }
      }
      gamma_mu0			<- gamma_hyp_prior$gamma_mu0
      gamma_K0			<- gamma_hyp_prior$gamma_K0
      gamma_nu			<- gamma_hyp_prior$gamma_nu
      gamma_V			  <- gamma_hyp_prior$gamma_V
    }
  }


  # Initialize emiss hyper prior
  if(missing(emiss_hyp_prior)){
    stop("The hyper-prior values for the Normal emission distribution(s) denoted by emiss_hyp_prior needs to be specified")
  }
  if (!is.mHMM_prior_emiss(emiss_hyp_prior) | !is.cont_f(emiss_hyp_prior)){
    stop("The input object specified for emiss_hyp_prior should be from the class mHMM_prior_emiss, obtained by using the function for continuous data with fixed emissions: prior_emiss_cont_f.")
  }
  if (emiss_hyp_prior$gen$m != m){
    stop("The number of states specified in m is not equal to the number of states specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cont_f.")
  }
  if (emiss_hyp_prior$gen$n_dep != n_dep){
    stop("The number of dependent variables specified in n_dep is not equal to the number of dependent variables specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cont_f.")
  }
  # emiss_mu0: a list containing n_dep matrices with in the first row the hypothesized mean values of the Normal emission
  # distributions in each of the states over the m columns. Subsequent rows contain the hypothesized regression
  # coefficients for covariates influencing the state dependent mean value of the normal distribution
  emiss_mu0	  <- emiss_hyp_prior$emiss_mu0
  emiss_V	    <- emiss_hyp_prior$emiss_V
  emiss_nu	  <- emiss_hyp_prior$emiss_nu
  emiss_K0    <- emiss_hyp_prior$emiss_K0



  # Define objects used to store data in mcmc algorithm, not returned ----------------------------
  # overall
  c <- llk <- numeric(1)
  sample_path <- lapply(n_vary, dif_matrix, cols = J)
  trans <- rep(list(vector("list", m)), n_subj)

  # gamma
  if(m > 1){
    gamma_int_mle_pooled  <- gamma_pooled_ll <- vector("list", m)
    gamma_c_int           <- rep(list(matrix(NA_real_, n_subj, (m-1))), m)
    gamma_mu_int_bar      <- gamma_V_int <- vector("list", m)
    gamma_naccept         <- matrix(0, n_subj, m)
    gamma_mu_prob_bar     <- rep(list(numeric(m)), m)
  } else if (m == 1){
    gamma_naccept         <- "With 1 state, this output object is obsolete"
  }

  # emiss
  cond_y <- lapply(rep(m, n_subj), nested_list, m = n_dep)
  emiss_c_mu <- rep(list(rep(list(matrix(NA_real_,ncol = 1, nrow = n_subj)),n_dep)), m)
  for(i in 1:m){
    for(q in 1:n_dep){
      emiss_c_mu[[i]][[q]][,1] <- start_val[[1 + q]][i,1]
    }
  }
  emiss_V_mu <- emiss_c_mu <- emiss_c_V <- rep(list(rep(list(NULL),n_dep)), m)
  ss_subj <- n_cond_y <- numeric(n_subj)


  # Define objects that are returned from mcmc algorithm ----------------------------
  # Define object for subject specific posterior density, put start values on first row
  if(length(start_val) != n_dep + 1 | depth(start_val) != 1){
    stop("The number of elements in the list start_val should be equal to 1 + the number of dependent variables,
         and should not contain nested lists (i.e., lists within lists)")
  }
  PD              <- list(trans_prob = matrix(NA_real_, nrow = J, ncol = m * m),
                          log_likl = matrix(NA_real_, nrow = J, ncol = 1))
  if(dim(start_val[[1]])[1] != m | dim(start_val[[1]])[2] != m){
    stop(paste0("Start values for the transition probability matrix contained in the first element of 'start_val' should be an m x m matrix, here ", m, " by ", m,"."))
  }
  PD$trans_prob[1, ] <- unlist(sapply(start_val, t))[1:(m*m)]
  if(m == 1){
    PD$trans_prob[,1] <- 1
  }
  colnames(PD$log_likl) <-  "LL"
  PD_subj				<- rep(list(PD), n_subj)

  # Define object for population posterior density (probabilities and regression coefficients parameterization )
  # gamma
  gamma_prob_bar		<- matrix(NA_real_, nrow = J, ncol = (m * m))
  colnames(gamma_prob_bar) <- paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = "")
  gamma_prob_bar[1,] <- PD$trans_prob[1, ]

  if(m > 1){
    gamma_int_bar				<- matrix(NA_real_, nrow = J, ncol = ((m-1) * m))
    colnames(gamma_int_bar) <- paste("int_S", rep(1:m, each = m-1), "toS", rep(2:m, m), sep = "")
    gamma_int_bar[1,] <- as.vector(t(prob_to_int(matrix(gamma_prob_bar[1,], byrow = TRUE, ncol = m, nrow = m))))
    gamma_V_idx <- which(paste("var_int_S", rep(1:m, each = (m-1)*(m-1)), "toS", rep(2:m, each=m-1), "_with_", "int_S", rep(1:m, each = (m-1)*(m-1)), "toS", rep(2:m, m), sep = "")
                         %in% paste0("var_int_S",rep(1:m,each=m-1),"toS",2:m,"_with_int_S",rep(1:m,each=m-1),"toS",2:m))
    gamma_V_int_bar <- matrix(NA_real_, nrow = J, ncol = ((m-1) * m))
    colnames(gamma_V_int_bar) <- paste0("var_int_S",rep(1:m,each=m-1),"toS",2:m)
    if(nx[1] > 1){
      gamma_cov_bar				<- matrix(NA_real_, nrow = J, ncol = ((m-1) * m) * (nx[1] - 1))
      colnames(gamma_cov_bar) <- paste( paste("cov", 1 : (nx[1] - 1), "_", sep = ""), "S", rep(1:m, each = (m-1) * (nx[1] - 1)), "toS", rep(2:m, m * (nx[1] - 1)), sep = "")
      gamma_cov_bar[1,] <- 0
    } else{
      gamma_cov_bar <- "No covariates where used to predict the transition probability matrix"
    }
    # Define object for subject specific posterior density (regression coefficients parameterization )
    gamma_int_subj			<- rep(list(gamma_int_bar), n_subj)
  } else if(m == 1){
    gamma_int_bar   <- "With 1 state, this output object is obsolete"
    gamma_V_int_bar <- "With 1 state, this output object is obsolete"
    gamma_cov_bar  <- "With 1 state, this output object is obsolete"
    gamma_int_subj  <- "With 1 state, this output object is obsolete"
  }


  # emiss
  emiss_mu			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("mu_", 1:m, sep = ""))))), n_dep)
  names(emiss_mu) <- dep_labels
  for(q in 1:n_dep){
    emiss_mu[[q]][1,] <- start_val[[q+1]][,1]
  }
  emiss_sd			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("sd_", 1:m, sep = ""))))), n_dep)
  names(emiss_sd) <- dep_labels
  for(q in 1:n_dep){
    emiss_sd[[q]][1,] <- start_val[[q+1]][,2]
  }


  # Put starting values in place for fist run forward algorithm
  emiss				<- start_val[2:(n_dep + 1)]
  emiss <- lapply(emiss, function(x) {
      x[,2] <- x[,2] ^2
      return(x)
  })
  if(sum(start_val[[1]] == 0) > 0){
    message("The starting vaules for your transition probability matrix gamma contains values equalling zero.
              This can cause estimation problems, please consider using (very small) non-zero values instead.")
  }
  gamma 			<- rep(start_val[1], n_subj)
  delta 			<- rep(list(solve(t(diag(m) - gamma[[1]] + 1), rep(1, m))), n_subj)



  # Start analysis --------------------------------------------
  # Run the MCMC algorithm
  itime <- proc.time()[3]
  if(show_progress == TRUE){
    cat("Progress of the Bayesian mHMM algorithm:", "\n")
    pb <- utils::txtProgressBar(min = 2, max = J, style = 3)
  }
  for (iter in 2 : J){

    if(m > 1){
      # For each subject, obtain sampled state sequence with subject individual parameters ----------
      for(s in 1:n_subj){
        # Run forward algorithm, obtain subject specific forward probabilities and log likelihood
        forward				<- cont_mult_fw_r_to_cpp(x = subj_data[[s]]$y, m = m, emiss = emiss, gamma = gamma[[s]], n_dep = n_dep, delta=NULL)
        alpha         <- abs(forward[[1]])
        if(sum(is.nan(alpha)) > 0){
          if(iter == 2){
            stop("The forward-backward algorithm ran into a fatal error during the first MCMC iteration. Most likely, starting values that do not match the data
               were specified, e.g., using means and (too small) standard deviations that do not sufficiently support the full range of the observed data
               (i.e., resulting in observation(s) that have an extremely small probability of being observed).")
          } else {
            "The forward-backward algorithm ran into a fatal error while running the MCMC algortihm. One possible cause could be the specification of hyper-prior
          distribution parameters that do not match the observed data. E.g., for continous data, setting prior possible values of the emission standard deviation
          by emiss_a0 and emiss_b0 too small. Please consider using differnt hyper-parameter values. "
          }
        }
        c             <- max(forward[[2]][, subj_data[[s]]$n_t])
        llk           <- c + log(sum(exp(forward[[2]][, subj_data[[s]]$n_t] - c)))
        PD_subj[[s]]$log_likl[iter, 1] <- llk

        # Using the forward probabilites, sample the state sequence in a backward manner.
        # In addition, saves state transitions in trans, and conditional observations within states in cond_y
        trans[[s]]					                  <- vector("list", m)
        sample_path[[s]][n_vary[[s]], iter] 	<- sample(1:m, 1, prob = c(alpha[, n_vary[[s]]]))
        for(t in (subj_data[[s]]$n_t - 1):1){
          sample_path[[s]][t,iter] 	              <- sample(1:m, 1, prob = (alpha[, t] * gamma[[s]][,sample_path[[s]][t + 1, iter]]))
          trans[[s]][[sample_path[[s]][t,iter]]]	<- c(trans[[s]][[sample_path[[s]][t, iter]]], sample_path[[s]][t + 1, iter])
        }
        for (i in 1:m){
          if(!is.null(trans[[s]][[i]])){
            trans[[s]][[i]] <- rev(trans[[s]][[i]])
          }
          for(q in 1:n_dep){
            cond_y[[s]][[q]][[i]] <- stats::na.omit(c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q]))
          }
        }
      }
    } else if(m == 1){
      for(s in 1:n_subj){
        allprobs <- all1(x = subj_data[[s]]$y, emiss = emiss[[s]], n_dep = n_dep, data_distr = data_distr)
        PD_subj[[s]]$log_likl[iter, 1] <- sum(log(stats::na.omit(allprobs)))
        for(q in 1:n_dep){
          cond_y[[s]][[q]][[i]] <- stats::na.omit(c(subj_data[[s]]$y[,q]))
        }
      }
    }
    # The remainder of the mcmc algorithm is state specific
    for(i in 1:m){

      # Obtain MLE of the covariance matrices and log likelihood of gamma and emiss at subject and population level -----------------
      # used to scale the proposal distribution of the RW Metropolis sampler

      # population level, transition matrix
      if(m > 1){
        trans_pooled			  <- factor(c(unlist(sapply(trans, "[[", i)), c(1:m)))
        gamma_mle_pooled		<- optim(gamma_int_mle0[i,], llmnl_int, Obs = trans_pooled,
                                   n_cat = m, method = "BFGS", hessian = FALSE,
                                   control = list(fnscale = -1))
        gamma_int_mle_pooled[[i]]  <- gamma_mle_pooled$par
        gamma_pooled_ll[[i]]			<- gamma_mle_pooled$value
      }

      # subject level
      for (s in 1:n_subj){
        wgt 				<- subj_data[[s]]$n_t / n_total

        # subject level, transition matrix
        if (m > 1){
          gamma_out					<- optim(gamma_int_mle_pooled[[i]], llmnl_int_frac, Obs = c(trans[[s]][[i]], c(1:m)),
                                 n_cat = m, pooled_likel = gamma_pooled_ll[[i]], w = gamma_w, wgt = wgt,
                                 method="BFGS", hessian = FALSE, control = list(fnscale = -1))
          if(gamma_out$convergence == 0){
            subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ]	<-
              mnlHess_int(int = gamma_out$par, Obs = c(trans[[s]][[i]], c(1:m)), n_cat =  m)
          } else {
            subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ]	<- diag(m-1)
          }
          # if this is first iteration, use MLE for current values RW metropolis sampler
          if (iter == 2){
            gamma_c_int[[i]][s,]		<- gamma_out$par
          }
        }
      }

      # Sample pouplaton values for gamma and conditional probabilities using Gibbs sampler -----------
      # gamma_mu0_n and gamma_mu_int_bar are matrices, with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
      if(m > 1){
        gamma_mu0_n           <- solve(t(xx[[1]]) %*% xx[[1]] + gamma_K0)  %*% (t(xx[[1]]) %*% gamma_c_int[[i]] + gamma_K0 %*% gamma_mu0[[i]])
        gamma_V_n             <- gamma_V + t(gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) %*% (gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) + t(gamma_mu0_n - gamma_mu0[[i]]) %*% gamma_K0 %*% (gamma_mu0_n - gamma_mu0[[i]])
        gamma_V_int[[i]]      <- solve(rwish(S = solve(gamma_V_n), v = gamma_nu + n_subj))
        gamma_mu_int_bar[[i]] <- gamma_mu0_n + solve(chol(t(xx[[1]]) %*% xx[[1]] + gamma_K0)) %*% matrix(rnorm((m - 1) * nx[1]), nrow = nx[1]) %*% t(solve(chol(solve(gamma_V_int[[i]]))))
        gamma_exp_int				  <- matrix(exp(c(0, gamma_mu_int_bar[[i]][1,] )), nrow  = 1)
        gamma_mu_prob_bar[[i]] 	<- gamma_exp_int / as.vector(gamma_exp_int %*% c(rep(1,(m))))
      }

      # sample fixed over subjects mean of the Normal emission distribution, and it's variance
      # note: since the emissions means are fixed over subjects, we can also not predict them using covariates!
      for(q in 1:n_dep){
        cond_y_pooled					         <- unlist(sapply(sapply(cond_y, "[[", q, simplify = FALSE), "[[", i))
        n_cond_y_pooled                <- length(cond_y_pooled)
        emiss_obs_mu                   <- mean(cond_y_pooled)
        emiss_mu0_n                    <- solve(n_cond_y_pooled + emiss_K0[[q]]) * (sum(cond_y_pooled) + emiss_K0[[q]] * emiss_mu0[[q]][,i])
        emiss_a_mu_n                   <- (emiss_nu[[q]] + n_cond_y_pooled) / 2
        emiss_b_mu_n                   <- (emiss_nu[[q]] * emiss_V[[q]][i]) / 2 +
          (matrix(cond_y_pooled - emiss_obs_mu, nrow = 1) %*% matrix(cond_y_pooled - emiss_obs_mu, ncol = 1) +
             solve(n_cond_y_pooled + emiss_K0[[q]]) * (n_cond_y_pooled * emiss_K0[[q]]) * (emiss_mu0[[q]][,i] - emiss_obs_mu)^2) / 2
        emiss_V_mu[[i]][[q]]       <- as.numeric(solve(stats::rgamma(1, shape = emiss_a_mu_n, rate = emiss_b_mu_n)))
        emiss_c_mu[[i]][[q]]	  <- emiss[[q]][i,1] <- emiss_mu0_n + rnorm(1, mean = 0, sd = sqrt(diag(emiss_V_mu[[i]][[q]] * solve(n_cond_y_pooled + emiss_K0[[q]]))))
        emiss_sd[[q]][iter, i] <- sqrt(emiss_V_mu[[i]][[q]])
        emiss[[q]][i,2] <- emiss_V_mu[[i]][[q]]
      }


      # Sample subject values for gamma and conditional probabilities using RW Metropolis sampler -----------
      if(m > 1){
        for (s in 1:n_subj){
          gamma_candcov_comb 			<- chol2inv(chol(subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ] + chol2inv(chol(gamma_V_int[[i]]))))
          gamma_RWout					    <- mnl_RW_once(int1 = gamma_c_int[[i]][s,], Obs = trans[[s]][[i]], n_cat = m, mu_int_bar1 = c(t(gamma_mu_int_bar[[i]]) %*% xx[[1]][s,]), V_int1 = gamma_V_int[[i]], scalar = gamma_scalar, candcov1 = gamma_candcov_comb)
          gamma[[s]][i,]  	      <- PD_subj[[s]]$trans_prob[iter, ((i-1) * m + 1) : ((i-1) * m + m)] <- (gamma_RWout$prob + .0001) / (1 + 0.0001 *m)
          gamma_naccept[s, i]			<- gamma_naccept[s, i] + gamma_RWout$accept
          gamma_c_int[[i]][s,]		<- gamma_RWout$draw_int
          gamma_int_subj[[s]][iter, (1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1))] <- gamma_c_int[[i]][s,]
          if(i == m){
            delta[[s]] 		<- solve(t(diag(m) - gamma[[s]] + 1), rep(1, m))
          }
        }
      }
    }


    # End of 1 MCMC iteration, save output values --------
    if(m > 1){
      gamma_int_bar[iter, ]				   	<- unlist(lapply(gamma_mu_int_bar, "[",1,))
      if(nx[1] > 1){
        gamma_cov_bar[iter, ]      	<- unlist(lapply(gamma_mu_int_bar, "[",-1,))
      }
      gamma_V_int_bar[iter, ] <- unlist(lapply(gamma_V_int, function(e) as.vector(t(e))))[gamma_V_idx]
      gamma_prob_bar[iter,]			<- unlist(gamma_mu_prob_bar)
    } else if(m == 1){
      gamma_prob_bar[iter,]			<- 1
    }

   for(q in 1:n_dep){
      emiss_mu[[q]][iter, ]	<- as.vector(unlist(lapply(
        lapply(emiss_c_mu, "[[", q), "[",1,)
      ))
    }

    if(show_progress == TRUE){
      utils::setTxtProgressBar(pb, iter)
    }
  }
  if(show_progress == TRUE){
    close(pb)
  }

  # End of function, return output values --------
  ctime = proc.time()[3]
  message(paste("Total time elapsed (hh:mm:ss):", hms(ctime-itime)))
  if(return_path == TRUE){
    out <- list(input = list(data_distr = data_distr, m = m, n_dep = n_dep, J = J,
                             burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                PD_subj = PD_subj, gamma_int_subj = gamma_int_subj,
                gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar,
                gamma_V_int_bar = gamma_V_int_bar,
                gamma_prob_bar = gamma_prob_bar,
                emiss_mu = emiss_mu, gamma_naccept = gamma_naccept,
                emiss_sd = emiss_sd,
                sample_path = sample_path)
  } else {
    out <- list(input = list(data_distr = data_distr, m = m, n_dep = n_dep, J = J,
                             burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                PD_subj = PD_subj, gamma_int_subj = gamma_int_subj,
                gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar,
                gamma_V_int_bar = gamma_V_int_bar,
                gamma_prob_bar = gamma_prob_bar,
                emiss_mu = emiss_mu, gamma_naccept = gamma_naccept,
                emiss_sd = emiss_sd)
  }
  class(out) <- append(class(out), c("mHMM_f", "cont_f"))
  return(out)
}
