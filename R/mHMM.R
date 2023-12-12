#' Multilevel hidden  Markov model using Bayesian estimation
#'
#' \code{mHMM} fits a multilevel (also known as mixed or random effects)
#' hidden Markov model (HMM) to intense longitudinal data with categorical
#' or continuous (i.e., normally distributed) observations of multiple
#' subjects using Bayesian estimation, and creates an object of class \code{mHMM}.
#' By using a multilevel framework, we allow for heterogeneity in the model
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
#'   columns contain the dependent variable(s). Note that in case of categorical
#'   dependent  variable(s), input are integers (i.e., whole numbers) that range
#'   from 1 to \code{q_emiss} (see below) and is numeric (i.e., not be a (set
#'   of) factor variable(s)). The total number of rows are equal to the sum over
#'   the number of observations of each subject, and the number of columns are
#'   equal to the number of dependent variables (\code{n_dep}) + 1. The number
#'   of observations can vary over subjects.
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
#'   number of dependent variables}
#'   \item{\code{q_emiss}: only to be specified if the data represents
#'   categorical data. Numeric vector with length \code{n_dep} denoting the
#'   number of observed categories for the categorical emission distribution
#'   for each of the dependent variables.}}
#' @param xx An optional list of (level 2) covariates to predict the transition
#'   matrix and/or the emission probabilities. Level 2 covariate(s) means that
#'   there is one observation per subject of each covariate. The first element
#'   in the list \code{xx} is used to predict the transition matrix. Subsequent
#'   elements in the list are used to predict the emission distribution of (each
#'   of) the dependent variable(s). Each element in the list is a matrix, with
#'   the number of rows equal to the number of subjects. The first column of
#'   each matrix represents the intercept, that is, a column only consisting of
#'   ones. Subsequent columns correspond to covariates used to predict the
#'   transition matrix / emission distribution. See \emph{Details} for more
#'   information on the use of covariates.
#'
#'   If \code{xx} is omitted completely, \code{xx} defaults to \code{NULL},
#'   resembling no covariates. Specific elements in the list can also be left
#'   empty (i.e., set to \code{NULL}) to signify that either the transition
#'   probability matrix or a specific emission distribution is not predicted by
#'   covariates.
#' @param start_val List containing the start values for the transition
#'   probability matrix gamma and the emission distribution(s). The first
#'   element of the list contains a \code{m} by \code{m} matrix with the start
#'   values for gamma. The subsequent \code{n_dep} elements each contain a
#'   matrix with the start values for the emission distribution(s). In case of
#'   categorical observations: a \code{m} by \code{q_emiss[k]} matrix denoting
#'   the probability of observing each categorical outcome (columns) within each
#'   state (rows). In case of continuous observations: a \code{m} by 2 matrix
#'   denoting the mean (first column) and standard deviation (second column) of
#'   the Normal emission distribution within each state (rows). In case of
#'   count observations: a \code{m} by 1 matrix denoting the mean (column) of
#'   the Poisson emission distribution within each state on the natural
#'   (i.e., Positive domain, not logarithmic) scale (rows). Note that
#'   \code{start_val} should not contain nested lists (i.e., lists within
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
#'   on the categorical or Normal emission distribution(s), generated by the
#'   function \code{\link{prior_emiss_cat}}, \code{\link{prior_emiss_cont}} or
#'  \code{\link{prior_emiss_count}}.
#'   Mandatory in case of continuous and count observations, optional in case of
#'   categorical observations.
#' @param gamma_sampler An optional object of the class \code{mHMM_pdRW_gamma}
#'   containing user specified settings for the proposal distribution of the
#'   random walk (RW) Metropolis sampler on the subject level transition
#'   probability matrix parameters, generated by the function
#'   \code{\link{pd_RW_gamma}}.
#' @param emiss_sampler An optional object of the class \code{mHMM_pdRW_emiss}
#'   containing user specified settings for the proposal distribution of the
#'   random walk (RW) Metropolis sampler on the subject level emission
#'   distribution(s) parameters, generated by the function
#'   \code{\link{pd_RW_emiss_cat}}. Only applicable in case of categorical
#'   observations.
#'
#' @return \code{mHMM} returns an object of class \code{mHMM}, which has
#'   \code{print} and \code{summary} methods to see the results.
#'   The object contains always the following components:
#'   \describe{
#'   \item{\code{PD_subj}}{A list containing one list per subject with the
#'   elements \code{trans_prob}, \code{cat_emiss} or \code{cont_emiss} in case
#'   of categorical or continuous observations, respectively, and
#'   \code{log_likl}, providing the subject parameter estimates over the
#'   iterations of the MCMC sampler. \code{trans_prob} relates to the transition
#'   probabilities gamma, \code{cat_emiss} to the categorical emission
#'   distribution (emission probabilities), \code{cont_emiss} to the continuous
#'   emission distributions (subsequently the the emission means and the (fixed
#'   over subjects) emission standard deviation), and \code{log_likl} to the log
#'   likelihood over the MCMC iterations. Iterations are contained in the rows,
#'   the parameters in the columns.}
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
#'   \item{\code{gamma_V_int_bar}}{A matrix containing the (co-)variance
#'   components for the subject-level intercepts
#'   of the multinomial logistic regression modeling the transition
#'   probabilities over the iterations of the hybrid Metropolis within Gibbs
#'   sampler. The iterations of the sampler are contained in the rows, and the
#'   columns contain the variance components for the subject level intercepts.}
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
#'   \item{\code{input}}{Overview of used input specifications: the distribution
#'   type of the observations \code{data_distr}, the number of
#'   states \code{m}, the number of used dependent variables \code{n_dep}, (in
#'   case of categorical observations) the number of output categories for each
#'   of the dependent variables \code{q_emiss}, the number of iterations
#'   \code{J} and the specified burn in period \code{burn_in} of the hybrid
#'   Metropolis within Gibbs sampler, the number of subjects \code{n_subj}, the
#'   observation length for each subject \code{n_vary}, and the column names of
#'   the dependent variables \code{dep_labels}.}
#'   \item{\code{sample_path}}{A list containing one matrix per subject with the
#'   sampled hidden state sequence over the hybrid Metropolis within Gibbs
#'   sampler. The time points of the dataset are contained in the rows, and the
#'   sampled paths over the iterations are contained in the columns. Only
#'   returned if \code{return_path = TRUE}. }
#' }
#'
#' Additionally, in case of categorical observations, the \code{mHMM} return object
#' contains:
#' \describe{
#'   \item{\code{emiss_prob_bar}}{A list containing one matrix per dependent
#'   variable, denoting the group level emission probabilities of each dependent
#'   variable over the iterations of the hybrid Metropolis within Gibbs sampler.
#'   The iterations of the sampler are contained in the rows of the matrix, and
#'   the columns contain the group level emission probabilities. If covariates
#'   were included in the analysis, the group level probabilities represent the
#'   predicted probability given that the covariate is at the average value for
#'   continuous covariates, or given that the covariate equals zero for
#'   dichotomous covariates.}
#'   \item{\code{emiss_int_bar}}{A list containing one matrix per dependent
#'   variable, denoting the group level intercepts of each dependent variable of
#'   the Multinomial logistic regression modeling the probabilities of the
#'   emission distribution over the iterations of the hybrid Metropolis within
#'   Gibbs sampler. The iterations of the sampler are contained in the rows of
#'   the matrix, and the columns contain the group level intercepts.}
#'   \item{\code{emiss_cov_bar}}{A list containing one matrix per dependent
#'   variable, denoting the group level regression coefficients of the
#'   Multinomial logistic regression predicting the emission probabilities within
#'   each of the dependent variables over the iterations of the hybrid
#'   Metropolis within Gibbs sampler. The iterations of the sampler are
#'   contained in the rows  of the matrix, and the columns contain the group
#'   level regression coefficients.}
#'   \item{\code{emiss_int_subj}}{A list containing one list per subject denoting
#'   the subject level intercepts of each dependent variable of the Multinomial
#'   logistic regression modeling the probabilities of the emission distribution
#'   over the iterations of the hybrid Metropolis within Gibbs sampler. Each
#'   lower level list contains one matrix per dependent variable, in which
#'   iterations of the sampler are contained in the rows, and the columns
#'   contain the subject level intercepts.}
#'   \item{\code{emiss_naccept}}{A list containing one matrix per dependent
#'   variable with the number of accepted draws at the subject level RW
#'   Metropolis step for each set of parameters of the emission distribution.
#'   The subjects are contained in the rows, and the columns of the matrix
#'   contain the sets of parameters.}}
#'
#' In case of continuous observations, the \code{mHMM} return object contains:
#' \describe{
#'   \item{\code{emiss_mu_bar}}{A list containing one matrix per dependent
#'   variable, denoting the group level means of the Normal emission
#'   distribution of each dependent variable over the iterations of the Gibbs
#'   sampler. The iterations of the sampler are contained in the rows of the
#'   matrix, and the columns contain the group level emission means. If
#'   covariates were included in the analysis, the group level means represent
#'   the predicted mean given that the covariate is at the average value for
#'   continuous covariates, or given that the covariate equals zero for
#'   dichotomous covariates.}
#'   \item{\code{emiss_varmu_bar}}{A list containing one matrix per dependent
#'   variable, denoting the variance between the subject level means of the
#'   Normal emission distributions. over the iterations of the Gibbs sampler. The
#'   iterations of the sampler are contained in the rows of the matrix, and the
#'   columns contain the group level variance in the mean.}
#'   \item{\code{emiss_sd_bar}}{A list containing one matrix per dependent
#'   variable, denoting the (fixed over subjects) standard deviation of the
#'   Normal emission distributions over the iterations of the Gibbs sampler. The
#'   iterations of the sampler are contained in the rows of the matrix, and the
#'   columns contain the group level emission variances.}
#'   \item{\code{emiss_cov_bar}}{A list containing one matrix per dependent
#'   variable, denoting the group level regression coefficients predicting the
#'   emission means within each of the dependent variables over the iterations
#'   of the Gibbs sampler. The iterations of the sampler are contained in the
#'   rows  of the matrix, and the columns contain the group level regression
#'   coefficients.}
#'   \item{\code{emiss_naccept}}{A list containing one matrix per dependent
#'   variable with the number of accepted draws at the subject level RW
#'   Metropolis step for each set of parameters of the emission distribution.
#'   The subjects are contained in the rows, and the columns of the matrix
#'   contain the sets of parameters.}
#'   \item{\code{input}}{Overview of used input specifications: the number of
#'   states \code{m}, the number of used dependent variables \code{n_dep}, the
#'   number of iterations \code{J} and the specified burn in period
#'   \code{burn_in} of the hybrid Metropolis within Gibbs sampler, the number of
#'   subjects \code{n_subj}, the observation length for each subject
#'   \code{n_vary}, and the column names of the dependent variables
#'   \code{dep_labels}.}}
#'
#' In case of count observations, the \code{mHMM} return object contains:
#' \describe{
#'   \item{\code{emiss_mu_bar}}{A list containing one matrix per dependent
#'   variable, denoting the group level means of the lognormal prior used
#'   in the emission distribution of each dependent variable over the
#'   MCMC iterations. The iterations of the sampler are contained in the rows
#'   of the matrix, and the columns contain the group-level lognormal emission
#'   logmeans in the logarithmic scale. If covariates were included in the
#'   analysis, the group level means represent the predicted mean given that the
#'   covariate is at the average value for continuous covariates, or given
#'   that the covariate equals zero for dichotomous covariates. Note that to
#'   obtain a corresponding median counts in the natual scale (i.e., Positive
#'   domain), the group-level emission logmeans have to be exponentiated).}
#'   \item{\code{emiss_varmu_bar}}{A list containing one matrix per dependent
#'   variable, denoting the logvariance between the subject-level means of the
#'   Poisson emission sampled from the lognormal prior distributions over the
#'   iterations of the MCMC sampler. The iterations of the sampler are
#'   contained in the rows of the matrix, and the columns contain the
#'   group-level logvariance between subject-specific means (in the
#'   logarithmic scale).}
#'   \item{\code{emiss_cov_bar}}{A list containing one matrix per dependent
#'   variable, denoting the group level regression coefficients predicting the
#'   emission means in the logarithmic scale within each of the dependent
#'   variables over the iterations of the MCMC sampler. The iterations of the
#'   sampler are contained in the rows  of the matrix, and the columns contain
#'   the group-level regression coefficients in the logarithmic scale.}
#'   \item{\code{input}}{Overview of used input specifications: the number of
#'   states \code{m}, the number of used dependent variables \code{n_dep}, the
#'   number of iterations \code{J} and the specified burn in period
#'   \code{burn_in} of the hybrid Metropolis within Gibbs sampler, the number of
#'   subjects \code{n_subj}, the observation length for each subject
#'   \code{n_vary}, and the column names of the dependent variables
#'   \code{dep_labels}.}}
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
#' ###### Example on package (categorical) example data, see ?nonverbal
#' \donttest{
#' # specifying general model properties:
#' m <- 2
#' n_dep <- 4
#' q_emiss <- c(3, 2, 3, 2)
#'
#' # specifying starting values
#' start_TM <- diag(.8, m)
#' start_TM[lower.tri(start_TM) | upper.tri(start_TM)] <- .2
#' start_EM <- list(matrix(c(0.05, 0.90, 0.05,
#'                           0.90, 0.05, 0.05), byrow = TRUE,
#'                         nrow = m, ncol = q_emiss[1]), # vocalizing patient
#'                  matrix(c(0.1, 0.9,
#'                           0.1, 0.9), byrow = TRUE, nrow = m,
#'                         ncol = q_emiss[2]), # looking patient
#'                  matrix(c(0.90, 0.05, 0.05,
#'                           0.05, 0.90, 0.05), byrow = TRUE,
#'                         nrow = m, ncol = q_emiss[3]), # vocalizing therapist
#'                  matrix(c(0.1, 0.9,
#'                           0.1, 0.9), byrow = TRUE, nrow = m,
#'                         ncol = q_emiss[4])) # looking therapist
#'
#' # Run a model without covariate(s):
#' # Note that for reasons of running time, J is set at a ridiculous low value.
#' # One would typically use a number of iterations J of at least 1000,
#' # and a burn_in of 200.
#' out_2st <- mHMM(s_data = nonverbal,
#'                 gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                 start_val = c(list(start_TM), start_EM),
#'                 mcmc = list(J = 11, burn_in = 5))
#'
#' out_2st
#' summary(out_2st)
#'
#' # plot the posterior densities for the transition and emission probabilities
#' plot(out_2st, component = "gamma", col =c("darkslategray3", "goldenrod"))
#'
#' # Run a model including a covariate (see ?nonverbal_cov) to predict the
#' # emission distribution for each of the 4 dependent variables:
#'
#' n_subj <- 10
#' xx_emiss <- rep(list(matrix(c(rep(1, n_subj),nonverbal_cov$std_CDI_change),
#'                             ncol = 2, nrow = n_subj)), n_dep)
#' xx <- c(list(matrix(1, ncol = 1, nrow = n_subj)), xx_emiss)
#' out_2st_c <- mHMM(s_data = nonverbal, xx = xx,
#'                   gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                   start_val = c(list(start_TM), start_EM),
#'                   mcmc = list(J = 11, burn_in = 5))
#'
#' }
#' ###### Example on categorical simulated data
#' # Simulate data for 10 subjects with each 100 observations:
#' n_t <- 100
#' n <- 10
#' m <- 2
#' n_dep <- 1
#' q_emiss <- 3
#' gamma <- matrix(c(0.8, 0.2,
#'                   0.3, 0.7), ncol = m, byrow = TRUE)
#' emiss_distr <- list(matrix(c(0.5, 0.5, 0.0,
#'                         0.1, 0.1, 0.8), nrow = m, ncol = q_emiss, byrow = TRUE))
#' data1 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                   gamma = gamma, emiss_distr = emiss_distr, var_gamma = .5, var_emiss = .5)
#'
#' # Specify remaining required analysis input (for the example, we use simulation
#' # input as starting values):
#' n_dep <- 1
#' q_emiss <- 3
#'
#' # Run the model on the simulated data:
#' out_2st_sim <- mHMM(s_data = data1$obs,
#'                  gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                  start_val = c(list(gamma), emiss_distr),
#'                  mcmc = list(J = 11, burn_in = 5))
#'
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
#' data_cont <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous', gen = list(m = m, n_dep = n_dep),
#'                   gamma = gamma, emiss_distr = emiss_distr, var_gamma = .1, var_emiss = c(.5, 0.01))
#'
#' # Specify hyper-prior for the continuous emission distribution
#' manual_prior_emiss <- prior_emiss_cont(
#'                         gen = list(m = m, n_dep = n_dep),
#'                         emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
#'                                          matrix(c(7, 8, 18), nrow = 1)),
#'                         emiss_K0 = list(1, 1),
#'                         emiss_V =  list(rep(100, m), rep(25, m)),
#'                         emiss_nu = list(1, 1),
#'                         emiss_a0 = list(rep(1, m), rep(1, m)),
#'                         emiss_b0 = list(rep(1, m), rep(1, m)))
#'
#' # Run the model on the simulated data:
#' # Note that for reasons of running time, J is set at a ridiculous low value.
#' # One would typically use a number of iterations J of at least 1000,
#' # and a burn_in of 200.
#' out_3st_cont_sim <- mHMM(s_data = data_cont$obs,
#'                          data_distr = 'continuous',
#'                          gen = list(m = m, n_dep = n_dep),
#'                          start_val = c(list(gamma), emiss_distr),
#'                          emiss_hyp_prior = manual_prior_emiss,
#'                          mcmc = list(J = 11, burn_in = 5))
#'
#' }
#'
#'
#' ###### Example on multivariate count data
#' # Simulate data with one covariate for each count dependent variable
#' \donttest{
#'
#' n_t     <- 200     # Number of observations on the dependent variable
#' m       <- 3        # Number of hidden states
#' n_dep   <- 3        # Number of dependent variables
#' n_subj  <- 30        # Number of subjects
#'
#' gamma   <- matrix(c(0.9, 0.05, 0.05,
#'                     0.2, 0.7, 0.1,
#'                     0.2,0.3, 0.5), ncol = m, byrow = TRUE)
#'
#' emiss_distr <- list(matrix(c(log(20),
#'                              log(10),
#'                              log(5)), nrow = m, byrow = TRUE),
#'                     matrix(c(log(15),
#'                              log(2),
#'                              log(5)), nrow = m, byrow = TRUE),
#'                     matrix(c(log(50),
#'                              log(3),
#'                              log(20)), nrow = m, byrow = TRUE))
#'
#' # Define list of vectors of covariate values
#' set.seed(42)
#' xx_vec <- c(list(NULL),rep(list(rnorm(n_subj,mean = 0, sd = 0.1)),3))
#'
#' # Define object beta with regression coefficients for the three dependent variables
#' beta      <- rep(list(NULL), n_dep+1)
#' beta[[2]] <- matrix(c(1,-1,0), byrow = TRUE, ncol = 1)
#' beta[[3]] <- matrix(c(2,0,-2), byrow = TRUE, ncol = 1)
#' beta[[4]] <- matrix(c(-1,3,1), byrow = TRUE, ncol = 1)
#'
#' # Create function to calculate between subject variance from logmu and logvar:
#' get_varmu <- function(lambda, logvar){
#'   logmu = log(lambda)
#'   abs(exp(logvar)-1)*exp(2*logmu+logvar)
#' }
#'
#' # Create function to calculate necessary logvar to get desired between
#' #   subject variance (depends also on the emission means chosen):
#' get_logvar <- function(mu, varmu){
#'   logmu = log(mu)
#'   log(0.5*exp(-2*logmu)*(exp(2*logmu) + sqrt(4*exp(2*logmu)*varmu+exp(4*logmu))))
#' }
#'
#' # Use the largest mean of each dependent variable:
#' logvar <- get_logvar(c(20,15,50), c(2,1.5,5)**2)
#'
#' # Simulate count data
#' data_count <- sim_mHMM(n_t = n_t, n = n_subj,
#'                        data_distr = "count", gen = list(m = m, n_dep = n_dep),
#'                        gamma = gamma, emiss_distr = emiss_distr, xx_vec = xx_vec, beta = beta,
#'                        var_gamma = 0.1, var_emiss = logvar, return_ind_par = TRUE)
#'
#' # Transition probabilities
#' start_gamma <- diag(0.8, m)
#' start_gamma[lower.tri(start_gamma) | upper.tri(start_gamma)] <- (1 - diag(start_gamma)) / (m - 1)
#'
#' # Emission distribution
#' start_emiss <- list(matrix(c(20,10, 5), nrow = m, byrow = TRUE),
#'                     matrix(c(15, 2, 5), nrow = m, byrow = TRUE),
#'                     matrix(c(50, 3,20), nrow = m, byrow = TRUE))
#'
#' # Specify list of matrices with covariates
#' set.seed(42)
#' xx_vec <- c(list(NULL),rep(list(matrix(cbind(1, rnorm(n_subj,mean = 0, sd = 0.1)), ncol = 2)),3))
#'
#' # Specify hyper-prior for the continuous emission distribution
#' manual_prior_emiss <- prior_emiss_count(
#'   gen = list(m = m, n_dep = n_dep),
#'   emiss_mu0 = list(matrix(log(c(20, 10, 5,
#'                                 1,  1, 1)), byrow = TRUE, ncol = m),
#'                    matrix(log(c(15, 2, 5,
#'                                 1, 1, 1)), byrow = TRUE, ncol = m),
#'                    matrix(log(c(50, 3, 20,
#'                                 1, 1,  1)), byrow = TRUE, ncol = m)),
#'   emiss_K0  = rep(list(rep(0.1, 2)),n_dep),
#'   emiss_nu  = rep(list(0.1),n_dep),
#'   emiss_V   = rep(list(rep(10, m)),n_dep),
#'   n_xx_emiss = rep(1,m)
#' )
#'
#' # Run model
#' # Note that for reasons of running time, J is set at a ridiculous low value.
#' # One would typically use a number of iterations J of at least 1000,
#' # and a burn_in of 200.
#' out_3st_count_sim <- mHMM(s_data = data_count$obs,
#'                           data_distr = 'count',
#'                           gen = list(m = m, n_dep = n_dep),
#'                           xx = xx_vec,
#'                           start_val = c(list(start_gamma), start_emiss),
#'                           emiss_hyp_prior = manual_prior_emiss,
#'                           mcmc = list(J = 11, burn_in = 5),
#'                           show_progress = TRUE)
#'}
#'
#' @export
#'
#'

mHMM <- function(s_data, data_distr = 'categorical', gen, xx = NULL, start_val, mcmc, return_path = FALSE, show_progress = TRUE,
                 gamma_hyp_prior = NULL, emiss_hyp_prior = NULL, gamma_sampler = NULL, emiss_sampler = NULL){
  # Initialize data -----------------------------------
  # dependent variable(s), sample size, dimensions gamma and conditional distribution
  if(sum(objects(gen) %in% "m") != 1 | sum(objects(gen) %in% "n_dep") != 1){
    stop("The input argument gen should contain the elements m and n_dep.")
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
  if(data_distr == 'categorical'){
    if( sum(objects(gen) %in% "q_emiss") != 1){
      stop("The input argument gen should contain the elements q_emiss when modelling categorical observations.")
    }
    q_emiss 		 <- gen$q_emiss
    if(length(q_emiss) != n_dep){
      stop("The lenght of q_emiss specifying the number of output categories for each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
    }
    if(sum(apply(as.matrix(s_data[,2:(n_dep + 1)], ncol = n_dep), 2, function(x) length(unique(x))) != q_emiss) > 0 |
       min(s_data[,2:(n_dep + 1)]) < 1 |
       sum(apply(as.matrix(s_data[,2:(n_dep + 1)], ncol = n_dep), 2, max) > q_emiss) > 0){
      stop("For categorical input variables, the values should be integers ranging from 1 to number of observed categories specified for each of the dependent variables in q_emiss")
    }
  } else {
    q_emiss <- rep(1, n_dep)
  }
  if (data_distr == 'categorical'){
    emiss_mhess   <- rep(list(NULL), n_dep)
    for(q in 1:n_dep){
      emiss_mhess[[q]] <- matrix(NA_real_, (q_emiss[q] - 1) * m, (q_emiss[q] - 1))
    }
  } else if (data_distr == 'count'){
    emiss_mhess     <- rep(list(rep(NA_real_,m)), n_dep)
  } else if (data_distr == 'continuous'){
    emiss_mhess     <- NULL
  }
  emiss_mle       <- rep(list(rep(NA_real_,m)), n_dep)
  emiss_converge  <- rep(list(rep(NA_real_,m)), n_dep)
  for(s in 1:n_subj){
    ypooled   <- rbind(ypooled, subj_data[[s]]$y)
    n_t       <- dim(subj_data[[s]]$y)[1]
    n_vary[s] <- n_t
    subj_data[[s]]	<- c(subj_data[[s]], n_t = n_t, list(gamma_mhess = matrix(NA_real_, (m - 1) * m, (m - 1)),
                                                        emiss_mle = emiss_mle,
                                                        emiss_mhess = emiss_mhess,
                                                        emiss_converge = emiss_converge))
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

  if(data_distr == 'categorical'){
    # Initialize emiss sampler
    if(is.null(emiss_sampler)){
      emiss_int_mle0 <- rep(list(NULL), n_dep)
      emiss_scalar	<- rep(list(NULL), n_dep)
      for(q in 1:n_dep){
        emiss_int_mle0[[q]] <- matrix(0, ncol = q_emiss[q] - 1, nrow = m)
        emiss_scalar[[q]] 	<- 2.93 / sqrt(q_emiss[q] - 1)
      }
      emiss_w		<- rep(list(.1), n_dep)
    } else {
      if (!is.mHMM_pdRW_emiss(emiss_sampler)){
        stop("The input object specified for emiss_sampler should be from the class mHMM_pdRW_emiss, obtained by using the function pd_RW_emiss_cat")
      }
      if (emiss_sampler$gen$m != m){
        stop("The number of states specified in m is not equal to the number of states specified when setting the proposal distribution of the RW Metropolis sampler on the emission distribution(s) using the function pd_RW_emiss_cat.")
      }
      if (emiss_sampler$gen$n_dep != n_dep){
        stop("The number of dependent variables specified in n_dep is not equal to the number of dependent variables specified when setting the proposal distribution of the RW Metropolis sampler on the emission distribution(s) using the function pd_RW_emiss_cat.")
      }
      if (sum(emiss_sampler$gen$q_emiss != q_emiss) > 0){
        stop("The number of number of observed categories for each of the dependent variable specified in q_emiss is not equal to q_emiss specified when setting the proposal distribution of the RW Metropolis sampler on the emission distribution(s) using the function pd_RW_emiss_cat.")
      }
      emiss_int_mle0	<- emiss_sampler$emiss_int_mle0
      emiss_scalar 	  <- emiss_sampler$emiss_scalar
      emiss_w    		  <- emiss_sampler$emiss_w
    }
  } else if(data_distr == 'count'){
    # Initialize emiss sampler
    if(is.null(emiss_sampler)) {
      emiss_scalar    <- rep(list(2.38), n_dep)
      emiss_w	        <- rep(list(.1), n_dep)
    } else {
      if (!is.mHMM_pdRW_emiss(emiss_sampler)){
        stop("The input object specified for emiss_sampler should be from the class mHMM_pdRW_emiss, obtained by using the function pd_RW_emiss_cat")
      }
      if (emiss_sampler$gen$m != m){
        stop("The number of states specified in m is not equal to the number of states specified when setting the proposal distribution of the RW Metropolis sampler on the emission distribution(s) using the function pd_RW_emiss_cat.")
      }
      if (emiss_sampler$gen$n_dep != n_dep){
        stop("The number of dependent variables specified in n_dep is not equal to the number of dependent variables specified when setting the proposal distribution of the RW Metropolis sampler on the emission distribution(s) using the function pd_RW_emiss_cat.")
      }
      emiss_scalar    <- emiss_sampler$emiss_scalar
      emiss_w         <- emiss_sampler$emiss_w
    }
  }


  # Initialize Gamma hyper prior
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

  # Initialize emiss hyper prior
  if(data_distr == 'categorical'){
    if(is.null(emiss_hyp_prior)){
      emiss_mu0	  <- rep(list(vector("list", m)), n_dep)
      emiss_nu	  <- rep(list(NULL), n_dep)
      emiss_V	    <- rep(list(NULL), n_dep)
      emiss_K0    <- rep(list(NULL), n_dep)
      for(q in 1:n_dep){
        for(i in 1:m){
          emiss_mu0[[q]][[i]]		<- matrix(0, ncol = q_emiss[q] - 1, nrow = nx[1 + q])
        }
        emiss_nu[[q]]		<- 3 + q_emiss[q] - 1
        emiss_V[[q]]		<- emiss_nu[[q]] * diag(q_emiss[q] - 1)
        emiss_K0[[q]]		<- diag(1, nx[1 + q])
      }
    } else {
      if (!is.mHMM_prior_emiss(emiss_hyp_prior) | !is.cat(emiss_hyp_prior)){
        stop("The input object specified for emiss_hyp_prior should be from the class mHMM_prior_emiss, obtained by using the function for categorical data: prior_emiss_cat.")
      }
      if (emiss_hyp_prior$gen$m != m){
        stop("The number of states specified in m is not equal to the number of states specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
      }
      if (emiss_hyp_prior$gen$n_dep != n_dep){
        stop("The number of dependent variables specified in n_dep is not equal to the number of dependent variables specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
      }
      if (sum(emiss_hyp_prior$gen$q_emiss != q_emiss) > 0){
        stop("The number of number of observed categories for each of the dependent variable specified in q_emiss is not equal to q_emiss specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
      }
      if(is.null(emiss_hyp_prior$n_xx_emiss) & sum(nx[2:n_dep1] > 1) > 0){
        stop("Covariates were specified to predict the emission distribution(s), but no covariates were specified when creating the informative hyper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
      }
      if(!is.null(emiss_hyp_prior$n_xx_emiss)){
        if(sum(emiss_hyp_prior$n_xx_emiss != (nx[2:n_dep1] - 1)) > 0){
          stop("The number of covariates specified to predict the emission distribution(s) is not equal to the number of covariates specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
        }
      }
      emiss_mu0	 <- emiss_hyp_prior$emiss_mu0
      emiss_nu   <- emiss_hyp_prior$emiss_nu
      emiss_V	   <- emiss_hyp_prior$emiss_V
      emiss_K0	 <- emiss_hyp_prior$emiss_K0
    }
  } else if(data_distr == 'continuous'){
    if(missing(emiss_hyp_prior)){
      stop("The hyper-prior values for the Normal emission distribution(s) denoted by emiss_hyp_prior needs to be specified")
    }
    if (!is.mHMM_prior_emiss(emiss_hyp_prior) | !is.cont(emiss_hyp_prior)){
      stop("The input object specified for emiss_hyp_prior should be from the class mHMM_prior_emiss, obtained by using the function for continuous data: prior_emiss_cont.")
    }
    if (emiss_hyp_prior$gen$m != m){
      stop("The number of states specified in m is not equal to the number of states specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cont.")
    }
    if (emiss_hyp_prior$gen$n_dep != n_dep){
      stop("The number of dependent variables specified in n_dep is not equal to the number of dependent variables specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cont.")
    }
    if(is.null(emiss_hyp_prior$n_xx_emiss) & sum(nx[2:n_dep1] > 1) > 0){
      stop("Covariates were specified to predict the emission distribution(s), but no covariates were specified when creating the informative hyper-prior distribution on the emission distribution(s) using the function prior_emiss_cont.")
    }
    if(!is.null(emiss_hyp_prior$n_xx_emiss)){
      if(sum(emiss_hyp_prior$n_xx_emiss != (nx[2:n_dep1] - 1)) > 0){
        stop("The number of covariates specified to predict the emission distribution(s) is not equal to the number of covariates specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
      }
    }
    # emiss_mu0: a list containing n_dep matrices with in the first row the hypothesized mean values of the Normal emission
    # distributions in each of the states over the m columns. Subsequent rows contain the hypothesized regression
    # coefficients for covariates influencing the state dependent mean value of the normal distribution
    emiss_mu0	  <- emiss_hyp_prior$emiss_mu0
    emiss_a0	  <- emiss_hyp_prior$emiss_a0
    emiss_b0	  <- emiss_hyp_prior$emiss_b0
    emiss_V	    <- emiss_hyp_prior$emiss_V
    emiss_nu	  <- emiss_hyp_prior$emiss_nu
    emiss_K0    <- emiss_hyp_prior$emiss_K0

  } else if(data_distr == 'count'){
    if(missing(emiss_hyp_prior)){
      stop("The hyper-prior values for the Normal emission distribution(s) denoted by emiss_hyp_prior needs to be specified")
    }
    if (!is.mHMM_prior_emiss(emiss_hyp_prior) | !is.count(emiss_hyp_prior)){
      stop("The input object specified for emiss_hyp_prior should be from the class mHMM_prior_emiss, obtained by using the function for continuous data: prior_emiss_cont.")
    }
    if (emiss_hyp_prior$gen$m != m){
      stop("The number of states specified in m is not equal to the number of states specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cont.")
    }
    if (emiss_hyp_prior$gen$n_dep != n_dep){
      stop("The number of dependent variables specified in n_dep is not equal to the number of dependent variables specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cont.")
    }
    if(is.null(emiss_hyp_prior$n_xx_emiss) & sum(nx[2:n_dep1] > 1) > 0){
      stop("Covariates were specified to predict the emission distribution(s), but no covariates were specified when creating the informative hyper-prior distribution on the emission distribution(s) using the function prior_emiss_cont.")
    }
    if(!is.null(emiss_hyp_prior$n_xx_emiss)){
      if(sum(emiss_hyp_prior$n_xx_emiss != (nx[2:n_dep1] - 1)) > 0){
        stop("The number of covariates specified to predict the emission distribution(s) is not equal to the number of covariates specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
      }
    }
    # emiss_mu0: a list containing n_dep matrices with in the first row the hypothesized mean values of the lognormal emission
    # distributions in each of the states over the m columns. Subsequent rows contain the hypothesized regression
    # coefficients for covariates influencing the state dependent mean value of the normal distribution. Note that
    # the hypothesized mean values are in the logarithmic scale.
    emiss_mu0	  <- emiss_hyp_prior$emiss_mu0
    emiss_nu	  <- emiss_hyp_prior$emiss_nu
    emiss_V	    <- emiss_hyp_prior$emiss_V
    emiss_K0	  <- emiss_hyp_prior$emiss_K0
  }



  # Define objects used to store data in mcmc algorithm, not returned ----------------------------
  # overall
  c <- llk <- numeric(1)
  sample_path <- lapply(n_vary, dif_matrix, cols = J)
  trans <- rep(list(vector("list", m)), n_subj)

  # gamma
  gamma_int_mle_pooled  <- gamma_pooled_ll <- vector("list", m)
  gamma_c_int           <- rep(list(matrix(NA_real_, n_subj, (m-1))), m)
  gamma_mu_int_bar      <- gamma_V_int <- vector("list", m)
  gamma_mu_prob_bar     <- rep(list(numeric(m)), m)
  gamma_naccept         <- matrix(0, n_subj, m)

  # emiss
  cond_y <- lapply(rep(m, n_subj), nested_list, m = n_dep)
  if(data_distr == 'categorical'){
    emiss_int_mle_pooled <- emiss_pooled_ll <- rep(list(vector("list", n_dep)), m)
    emiss_c_int <- rep(list(lapply(q_emiss - 1, dif_matrix, rows = n_subj)), m)
    emiss_mu_int_bar <- emiss_V_int <- rep(list(vector("list", n_dep)), m)
    emiss_mu_prob_bar <- rep(list(lapply(q_emiss, dif_vector)), m)
    emiss_naccept <- rep(list(matrix(0, n_subj, m)), n_dep)
  } else if(data_distr == 'continuous'){
    emiss_c_mu <- rep(list(rep(list(matrix(NA_real_,ncol = 1, nrow = n_subj)),n_dep)), m)
    for(i in 1:m){
      for(q in 1:n_dep){
        emiss_c_mu[[i]][[q]][,1] <- start_val[[1 + q]][i,1]
      }
    }
    emiss_V_mu <- emiss_c_mu_bar <- emiss_c_V <- rep(list(rep(list(NULL),n_dep)), m)
    ss_subj <- n_cond_y <- numeric(n_subj)
  } else if(data_distr == 'count'){
    cond_y_pooled <- rep(list(rep(list(NULL),n_dep)), m)
    emiss_c_mu <- rep(list(rep(list(matrix(NA_real_,ncol = 1, nrow = n_subj)),n_dep)), m)
    for(i in 1:m){
      for(q in 1:n_dep){
        emiss_c_mu[[i]][[q]][,1] <- start_val[[1 + q]][i,1]
      }
    }
    emiss_V_mu <- emiss_c_mu_bar <- emiss_c_V <- rep(list(rep(list(NULL),n_dep)), m)
    emiss_mle_pooled <- emiss_pooled_ll <- rep(list(rep(list(NULL),n_dep)), m)
    ss_subj <- n_cond_y <- numeric(n_subj)
    emiss_naccept <- rep(list(matrix(0, n_subj, m)), n_dep)
  }


  # Define objects that are returned from mcmc algorithm ----------------------------
  # Define object for subject specific posterior density, put start values on first row
  if(length(start_val) != n_dep + 1 | depth(start_val) != 1){
    stop("The number of elements in the list start_val should be equal to 1 + the number of dependent variables,
         and should not contain nested lists (i.e., lists within lists)")
  }
  PD              <- list(trans_prob = matrix(NA_real_, nrow = J, ncol = m * m),
                          cat_emiss = if(data_distr == 'categorical'){
                            matrix(NA_real_, nrow = J, ncol = sum(m * q_emiss))} else {
                              NULL
                            },
                          cont_emiss = if(data_distr == 'continuous') {
                            matrix(NA_real_, nrow = J, ncol = (n_dep * 2 * m))} else {
                              NULL
                            },
                          count_emiss = if(data_distr == 'count') {
                            matrix(NA_real_, nrow = J, ncol = (n_dep * m))} else {
                              NULL
                            },
                          log_likl = matrix(NA_real_, nrow = J, ncol = 1))
  if(dim(start_val[[1]])[1] != m | dim(start_val[[1]])[2] != m){
    stop(paste0("Start values for the transition probability matrix contained in the first element of 'start_val' should be an m x m matrix, here ", m, " by ", m,"."))
  }
  PD$trans_prob[1, ] <- unlist(sapply(start_val, t))[1:(m*m)]
  if (data_distr == 'categorical'){
    if(sum(sapply(start_val, dim)[1, 2:(n_dep+1)] != rep(m, n_dep)) > 0 |
       sum(sapply(start_val, dim)[2, 2:(n_dep+1)] != q_emiss) > 0){
      stop("Start values for the emission distribuitons for categorical observations contained in 'start_val' should be one matrix per dependent variable, each with m rows and q_emiss[k] columns")
    }
    PD_cat_emiss_names   <- character()
    for(q in 1:n_dep){
      PD_cat_emiss_names <- c(PD_cat_emiss_names, paste("dep", q, "_S", rep(1:m, each = q_emiss[q]), "_emiss", rep(1:q_emiss[q], m), sep = ""))
    }
    colnames(PD$cat_emiss)    <- PD_cat_emiss_names
    for(q in 1:n_dep){
      PD$cat_emiss[1, (sum(c(0,q_emiss)[1 : q] * m) + 1):sum(q_emiss[1 : q] * m)] <- as.vector(t(start_val[[q + 1]]))
    }
  } else if(data_distr == 'continuous'){
    if(sum(sapply(start_val, dim)[1, 2:(n_dep+1)] != rep(m, n_dep)) > 0 |
       sum(sapply(start_val, dim)[2, 2:(n_dep+1)] != rep(2, n_dep)) > 0){
      stop("Start values for the emission distribuitons for continuous observations contained in 'start_val' should be one matrix per dependent variable, each with m rows and 2 columns")
    }
    colnames(PD$cont_emiss) <- c(paste("dep", rep(1:n_dep, each = m), "_S", rep(1:m), "_mu", sep = ""),
                                 paste("dep", rep(1:n_dep, each = m), "_S", rep(1:m), "_sd", sep = ""))
    for(q in 1:n_dep){
      PD$cont_emiss[1, ((q-1) * m + 1):(q * m)] <- start_val[[q + 1]][,1]
      PD$cont_emiss[1,  (n_dep * m + (q-1) * m + 1):(n_dep * m + q * m)] <- start_val[[q + 1]][,2]^2
    }
  } else if(data_distr == 'count'){
    if(sum(sapply(start_val, dim)[1, 2:(n_dep+1)] != rep(m, n_dep)) > 0){
      stop("Start values for the emission distribuitons for continuous observations contained in 'start_val' should be one matrix per dependent variable, each with m rows and 1 column")
    }
    colnames(PD$count_emiss) <- c(paste("dep", rep(1:n_dep, each = m), "_S", rep(1:m), "_mu", sep = ""))
    for(q in 1:n_dep){
      PD$count_emiss[1, ((q-1) * m + 1):(q * m)] <- start_val[[q + 1]][,1]
    }
  }
  colnames(PD$log_likl) <-  "LL"
  PD_subj				<- rep(list(PD), n_subj)

  # Define object for population posterior density (probabilities and regression coefficients parameterization )
  # gamma
  gamma_prob_bar		<- matrix(NA_real_, nrow = J, ncol = (m * m))
  colnames(gamma_prob_bar) <- paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = "")
  gamma_prob_bar[1,] <- PD$trans_prob[1, ]
  gamma_int_bar				<- matrix(NA_real_, nrow = J, ncol = ((m-1) * m))
  colnames(gamma_int_bar) <- paste("int_S", rep(1:m, each = m-1), "toS", rep(2:m, m), sep = "")
  gamma_int_bar[1,] <- as.vector(t(prob_to_int(matrix(gamma_prob_bar[1,], byrow = TRUE, ncol = m, nrow = m))))
  gamma_V_int_bar <- matrix(NA_real_, nrow = J, ncol = ((m-1) * (m-1) * m))
  colnames(gamma_V_int_bar) <- paste("var_int_S", rep(1:m, each = (m-1)*(m-1)), "toS", rep(2:m, each=m-1), "_with_", "int_S", rep(1:m, each = (m-1)*(m-1)), "toS", rep(2:m, m), sep = "")
  gamma_V_int_bar[1,] <- unlist(lapply(gamma_V, function(e) as.vector(t(e))))
  if(nx[1] > 1){
    gamma_cov_bar				<- matrix(NA_real_, nrow = J, ncol = ((m-1) * m) * (nx[1] - 1))
    colnames(gamma_cov_bar) <- paste( paste("cov", 1 : (nx[1] - 1), "_", sep = ""), "S", rep(1:m, each = (m-1) * (nx[1] - 1)), "toS", rep(2:m, m * (nx[1] - 1)), sep = "")
    gamma_cov_bar[1,] <- 0
  } else{
    gamma_cov_bar <- "No covariates where used to predict the transition probability matrix"
  }
  # Define object for subject specific posterior density (regression coefficients parameterization )
  gamma_int_subj			<- rep(list(gamma_int_bar), n_subj)

  # emiss
  if(data_distr == 'categorical'){
    emiss_prob_bar			<- lapply(q_emiss * m, dif_matrix, rows = J)
    names(emiss_prob_bar) <- dep_labels
    for(q in 1:n_dep){
      colnames(emiss_prob_bar[[q]]) <- paste("Emiss", rep(1:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q]), sep = "")
      start <- c(0, q_emiss * m)
      emiss_prob_bar[[q]][1,] <- PD$cat_emiss[1, (sum(c(0,q_emiss)[1 : q] * m) + 1):sum(q_emiss[1 : q] * m)]
    }
    emiss_int_bar			<- lapply((q_emiss-1) * m, dif_matrix, rows = J)
    names(emiss_int_bar) <- dep_labels
    for(q in 1:n_dep){
      colnames(emiss_int_bar[[q]]) <-  paste("int_Emiss", rep(2:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q] - 1), sep = "")
      emiss_int_bar[[q]][1,] <- as.vector(t(prob_to_int(matrix(emiss_prob_bar[[q]][1,], byrow = TRUE, ncol = q_emiss[q], nrow = m))))
    }
    if(sum(nx[-1]) > n_dep){
      emiss_cov_bar			<- lapply((q_emiss-1) * m * (nx[-1] - 1 ), dif_matrix, rows = J)
      names(emiss_cov_bar) <- dep_labels
      for(q in 1:n_dep){
        if(nx[1 + q] > 1){
          colnames(emiss_cov_bar[[q]]) <-  paste( paste("cov", 1 : (nx[1 + q] - 1), "_", sep = ""), "emiss", rep(2:q_emiss[q], m * (nx[1 + q] - 1)), "_S", rep(1:m, each = (q_emiss[q] - 1) * (nx[1 + q] - 1)), sep = "")
          emiss_cov_bar[[q]][1,] <- 0
        } else {
          emiss_cov_bar[[q]] <- "No covariates where used to predict the emission probabilities for this outcome"
        }
      }
    } else{
      emiss_cov_bar <- "No covariates where used to predict the emission probabilities"
    }
    # Define object for subject specific posterior density (regression coefficients parameterization )
    emiss_int_subj			<- rep(list(emiss_int_bar), n_subj)
  } else if(data_distr == 'continuous'){
    emiss_mu_bar			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("mu_", 1:m, sep = ""))))), n_dep)
    names(emiss_mu_bar) <- dep_labels
    for(q in 1:n_dep){
      emiss_mu_bar[[q]][1,] <- PD$cont_emiss[1, ((q-1) * m + 1):(q * m)]
    }
    emiss_varmu_bar			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("varmu_", 1:m, sep = ""))))), n_dep)
    names(emiss_varmu_bar) <- dep_labels
    emiss_sd_bar			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("sd_", 1:m, sep = ""))))), n_dep)
    names(emiss_sd_bar) <- dep_labels
    for(q in 1:n_dep){
      emiss_sd_bar[[q]][1,] <- PD$cont_emiss[1,  (n_dep * m + (q-1) * m + 1):(n_dep * m + q * m)]
    }
    if(sum(nx[-1]) > n_dep){
      emiss_cov_bar			<- lapply(m * (nx[-1] - 1 ), dif_matrix, rows = J)
      names(emiss_cov_bar) <- dep_labels
      for(q in 1:n_dep){
        if(nx[1 + q] > 1){
          colnames(emiss_cov_bar[[q]]) <-  paste( paste("cov", 1 : (nx[1+q] - 1), "_", sep = ""), "mu_S", rep(1:m, each = (nx[1 + q] - 1)), sep = "")
          emiss_cov_bar[[q]][1,] <- 0
        } else {
          emiss_cov_bar[[q]] <- "No covariates where used to predict the emission probabilities for this outcome"
        }
      }
    } else{
      emiss_cov_bar <- "No covariates where used to predict the emission probabilities"
    }
  } else if (data_distr == 'count'){
    emiss_mu_bar			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("mu_", 1:m, sep = ""))))), n_dep)
    names(emiss_mu_bar) <- dep_labels
    for(q in 1:n_dep){
      emiss_mu_bar[[q]][1,] <- log(start_val[[q + 1]][,1])
    }
    emiss_varmu_bar			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("varmu_", 1:m, sep = ""))))), n_dep)
    names(emiss_varmu_bar) <- dep_labels
    if(sum(nx[-1]) > n_dep){
      emiss_cov_bar			<- lapply(m * (nx[-1] - 1 ), dif_matrix, rows = J)
      names(emiss_cov_bar) <- dep_labels
      for(q in 1:n_dep){
        if(nx[1 + q] > 1){
          colnames(emiss_cov_bar[[q]]) <-  paste( paste("cov", 1 : (nx[1+q] - 1), "_", sep = ""), "mu_S", rep(1:m, each = (nx[1 + q] - 1)), sep = "")
          emiss_cov_bar[[q]][1,] <- 0
        } else {
          emiss_cov_bar[[q]] <- "No covariates where used to predict the emission probabilities for this outcome"
        }
      }
    } else{
      emiss_cov_bar <- "No covariates where used to predict the emission probabilities"
    }
  }

  # Put starting values in place for fist run forward algorithm
  if(data_distr == 'categorical'){
    if(sum(unlist(start_val[2:(n_dep + 1)]) == 0) > 0){
      message("The starting vaules for your categorical emission distribution(s) contain values equalling zero.
              This can cause estimation problems, please consider using (very small) non-zero values instead.")
    }
  }
  emiss				<- rep(list(start_val[2:(n_dep + 1)]), n_subj)
  if(data_distr == 'continuous'){
    emiss <- lapply(emiss, function(x) {
      lapply(x, function(x) {
        x[,2] <- x[,2] ^2
        return(x)
      })
    })
  }
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

    # For each subject, obtain sampled state sequence with subject individual parameters ----------
    for(s in 1:n_subj){
      # Run forward algorithm, obtain subject specific forward probabilities and log likelihood
      if(data_distr == 'categorical'){
        forward				<- cat_mult_fw_r_to_cpp(x = subj_data[[s]]$y, m = m, emiss = emiss[[s]], gamma = gamma[[s]], n_dep = n_dep, delta=NULL)
      } else if(data_distr == 'continuous'){
        forward				<- cont_mult_fw_r_to_cpp(x = subj_data[[s]]$y, m = m, emiss = emiss[[s]], gamma = gamma[[s]], n_dep = n_dep, delta=NULL)
      } else if(data_distr == 'count'){
        forward				<- count_mult_fw_r_to_cpp(x = subj_data[[s]]$y, m = m, emiss = emiss[[s]], gamma = gamma[[s]], n_dep = n_dep, delta=NULL)
      }
      alpha         <- abs(forward[[1]])
      if(sum(is.nan(alpha)) > 0){
        if(iter == 2){
          if(data_distr == 'continuous'){
            stop("The forward-backward algorithm ran into a fatal error during the first MCMC iteration. Most likely, starting values that do not match the data
                 were specified, e.g., using means and (too small) standard deviations that do not sufficiently support the full range of the observed data
                 (i.e., resulting in observation(s) that have an extremely small probability of being observed).")
          } else if (data_distr == 'count') {
            stop("The forward-backward algorithm ran into a fatal error during the first MCMC iteration. Most likely, starting values that do not match the data
                 were specified, e.g., using poisson means that do not sufficiently support the full range of the observed data
                 (i.e., resulting in observation(s) that have an extremely small probability of being observed). Notice that the starting valyes for the Poisson
                 parameters have to be specified in the Natural domain (i.e., not in the logarithmic domain).")
          } else if(data_distr == 'categorical'){
            stop("The forward-backward algorithm ran into a fatal error during the first MCMC iteration. Most likely, starting values that do not match the data
                 were specified, e.g., using emission categorical probabilities that do no sufficently support the full range of the observed data
                 (i.e., resulting in observation(s) that have an extremely small probability of being observed).")
          }
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

    # The remainder of the mcmc algorithm is state specific
    for(i in 1:m){

      # Obtain MLE of the covariance matrices and log likelihood of gamma and emiss at subject and population level -----------------
      # used to scale the proposal distribution of the RW Metropolis sampler

      # population level, transition matrix
      trans_pooled			  <- factor(c(unlist(sapply(trans, "[[", i)), c(1:m)))
      gamma_mle_pooled		<- optim(gamma_int_mle0[i,], llmnl_int, Obs = trans_pooled,
                                   n_cat = m, method = "BFGS", hessian = FALSE,
                                   control = list(fnscale = -1))
      gamma_int_mle_pooled[[i]]  <- gamma_mle_pooled$par
      gamma_pooled_ll[[i]]			<- gamma_mle_pooled$value

      # population level, conditional probabilities, separate for each dependent variable
      if(data_distr == 'categorical'){
        for(q in 1:n_dep){
          cond_y_pooled					     <- c(unlist(sapply(sapply(cond_y, "[[", q, simplify = FALSE), "[[", i)), c(1:q_emiss[q]))
          emiss_mle_pooled		       <- optim(emiss_int_mle0[[q]][i,], llmnl_int, Obs = cond_y_pooled,
                                              n_cat = q_emiss[q], method = "BFGS", hessian = FALSE,
                                              control = list(fnscale = -1))
          emiss_int_mle_pooled[[i]][[q]]  <- emiss_mle_pooled$par
          emiss_pooled_ll[[i]][[q]]				<- emiss_mle_pooled$value
        }
      } else if (data_distr == 'count') {
        for(q in 1:n_dep){
          cond_y_pooled					      <- numeric()
          for(s in 1:n_subj){
            cond_y_pooled             <- c(cond_y_pooled, cond_y[[s]][[q]][[i]])
          }
          if(iter == 2) {
            emiss_c_mu_bar[[i]][[q]] <- log(start_val[[q+1]][i,1])
          }
          emiss_mle_pooled[[i]][[q]]  <- mean(c(cond_y_pooled, round(exp(emiss_c_mu_bar[[i]][[q]][1]),0)), na.rm=TRUE)
          if(emiss_mle_pooled[[i]][[q]] == 0){
            emiss_mle_pooled[[i]][[q]] <- 0.01
          }
          emiss_pooled_ll[[i]][[q]]	<- llpois(lambda = emiss_mle_pooled[[i]][[q]], Obs = c(cond_y_pooled, round(exp(emiss_c_mu_bar[[i]][[q]][1]),0)))
        }
      }

      # subject level
      for (s in 1:n_subj){
        wgt 				<- subj_data[[s]]$n_t / n_total

        # subject level, transition matrix
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

        # subject level, conditional probabilities, separate for each dependent variable
        if(data_distr == 'categorical'){
          for(q in 1:n_dep){
            emiss_out				<- optim(emiss_int_mle_pooled[[i]][[q]], llmnl_int_frac, Obs = c(cond_y[[s]][[q]][[i]], c(1:q_emiss[q])),
                                  n_cat = q_emiss[q], pooled_likel = emiss_pooled_ll[[i]][[q]],
                                  w = emiss_w[[q]], wgt = wgt, method = "BFGS", hessian = FALSE, control = list(fnscale = -1))
            if(emiss_out$convergence == 0){
              subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ]		<-
                mnlHess_int(int = emiss_out$par, Obs = c(cond_y[[s]][[q]][[i]], c(1:q_emiss[q])), n_cat =  q_emiss[q])
            } else {
              subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ]	<- diag(q_emiss[q] - 1)
            }
            # if this is first iteration, use MLE for current values RW metropolis sampler
            if (iter == 2){
              emiss_c_int[[i]][[q]][s,]	<- emiss_out$par
            }
          }
        } else if (data_distr == 'count'){
          for(q in 1:n_dep){
            emiss_out	<- optim(log(emiss_mle_pooled[[i]][[q]]), llpois_frac_log, Obs = c(cond_y[[s]][[q]][[i]], round(exp(emiss_c_mu_bar[[i]][[q]][1]),0)),
                               pooled_likel = emiss_pooled_ll[[i]][[q]],
                               w = emiss_w[[q]], wgt = wgt,
                               method = "BFGS",
                               hessian = TRUE,
                               control = list(fnscale = -1))
            if(emiss_out$convergence == 0){
              subj_data[[s]]$emiss_converge[[q]][i]   <- 1
              subj_data[[s]]$emiss_mle[[q]][i]        <- exp(emiss_out$par)
              subj_data[[s]]$emiss_mhess[[q]][i]      <- as.numeric(-emiss_out$hessian)
            } else {
              subj_data[[s]]$emiss_converge[[q]][i]   <- 0
              subj_data[[s]]$emiss_mle[[q]][i]        <- 0
              subj_data[[s]]$emiss_mhess[[q]][i]	    <- 0 # check this
            }
          }
        }
      }

      # Sample pouplaton values for gamma and conditional probabilities using Gibbs sampler -----------
      # gamma_mu0_n and gamma_mu_int_bar are matrices, with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
      gamma_mu0_n           <- solve(t(xx[[1]]) %*% xx[[1]] + gamma_K0)  %*% (t(xx[[1]]) %*% gamma_c_int[[i]] + gamma_K0 %*% gamma_mu0[[i]])
      gamma_V_n             <- gamma_V + t(gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) %*% (gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) + t(gamma_mu0_n - gamma_mu0[[i]]) %*% gamma_K0 %*% (gamma_mu0_n - gamma_mu0[[i]])
      gamma_V_int[[i]]      <- solve(rwish(S = solve(gamma_V_n), v = gamma_nu + n_subj))
      gamma_mu_int_bar[[i]] <- gamma_mu0_n + solve(chol(t(xx[[1]]) %*% xx[[1]] + gamma_K0)) %*% matrix(rnorm((m - 1) * nx[1]), nrow = nx[1]) %*% t(solve(chol(solve(gamma_V_int[[i]]))))
      gamma_exp_int				  <- matrix(exp(c(0, gamma_mu_int_bar[[i]][1,] )), nrow  = 1)
      gamma_mu_prob_bar[[i]] 	<- gamma_exp_int / as.vector(gamma_exp_int %*% c(rep(1,(m))))

      if(data_distr == 'categorical'){
        for(q in 1:n_dep){
          emiss_mu0_n                 <- solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% (t(xx[[1 + q]]) %*% emiss_c_int[[i]][[q]] + emiss_K0[[q]] %*% emiss_mu0[[q]][[i]])
          emiss_V_n                   <- emiss_V[[q]] + t(emiss_c_int[[i]][[q]] - xx[[1 + q]] %*% emiss_mu0_n) %*% (emiss_c_int[[i]][[q]] - xx[[1 + q]] %*% emiss_mu0_n) + t(emiss_mu0_n - emiss_mu0[[q]][[i]]) %*% emiss_K0[[q]] %*% (emiss_mu0_n - emiss_mu0[[q]][[i]])
          emiss_V_int[[i]][[q]]       <- solve(rwish(S = solve(emiss_V_n), v = emiss_nu[[q]] + n_subj))
          emiss_mu_int_bar[[i]][[q]]	 <- emiss_mu0_n + solve(chol(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]])) %*% matrix(rnorm((q_emiss[q] - 1) * nx[1 + q]), nrow = nx[1 + q]) %*% t(solve(chol(solve(emiss_V_int[[i]][[q]]))))
          emiss_exp_int				       <- matrix(exp(c(0, emiss_mu_int_bar[[i]][[q]][1, ])), nrow  = 1)
          emiss_mu_prob_bar[[i]][[q]] <- emiss_exp_int / as.vector(emiss_exp_int %*% c(rep(1, (q_emiss[q]))))
        }
      } else if(data_distr == 'continuous'){
        # sample population mean (and regression parameters if covariates) of the Normal emission distribution, and it's variance (so the variance between the subject specific means)
        # note: the posterior is thus one of a Bayesian linear regression because of the optional regression parameters
        for(q in 1:n_dep){
          emiss_mu0_n                    <- solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% (t(xx[[1 + q]]) %*% emiss_c_mu[[i]][[q]] + emiss_K0[[q]] %*% emiss_mu0[[q]][,i])
          emiss_a_mu_n                   <- (emiss_K0[[q]] + n_subj) / 2
          emiss_b_mu_n                   <- (emiss_nu[[q]] * emiss_V[[q]][i]) / 2 + (t(emiss_c_mu[[i]][[q]]) %*% emiss_c_mu[[i]][[q]] +
                                                                                       t(emiss_mu0[[q]][,i]) %*% emiss_K0[[q]] %*% emiss_mu0[[q]][,i] -
                                                                                       t(emiss_mu0_n) %*% (t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% emiss_mu0_n) / 2
          emiss_V_mu[[i]][[q]]       <- as.numeric(solve(stats::rgamma(1, shape = emiss_a_mu_n, rate = emiss_b_mu_n)))
          emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(diag(emiss_V_mu[[i]][[q]] * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]))))
        }
      } else if(data_distr == 'count'){
        # sample population mean (and regression parameters if covariates) of the lognormal emission distribution, and it's variance (so the variance between the subject specific means)
        # note: the posterior is thus one of a Bayesian linear regression because of the optional regression parameters
        for(q in 1:n_dep){
          emiss_mu0_n                    <- solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% (t(xx[[1 + q]]) %*% log(emiss_c_mu[[i]][[q]]) + emiss_K0[[q]] %*% emiss_mu0[[q]][,i])
          emiss_a_mu_n                   <- (emiss_K0[[q]] + n_subj) / 2
          emiss_b_mu_n                   <- (emiss_nu[[q]] * emiss_V[[q]][i]) / 2 + (t(log(emiss_c_mu[[i]][[q]])) %*% log(emiss_c_mu[[i]][[q]]) +
                                                                                       t(emiss_mu0[[q]][,i]) %*% emiss_K0[[q]] %*% emiss_mu0[[q]][,i] -
                                                                                       t(emiss_mu0_n) %*% (t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% emiss_mu0_n) / 2
          emiss_V_mu[[i]][[q]]       <- as.numeric(solve(rgamma(1, shape = emiss_a_mu_n, rate = emiss_b_mu_n)))
          emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(diag(emiss_V_mu[[i]][[q]] * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]])))) # check
        }
      }



      # Sample subject values for gamma and conditional probabilities using RW Metropolis sampler -----------
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

      if (data_distr == 'categorical'){
        for (s in 1:n_subj){
          start <- c(0, q_emiss * m)
          for(q in 1:n_dep){
            emiss_candcov_comb		     <- chol2inv(chol(subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ] + chol2inv(chol(emiss_V_int[[i]][[q]]))))
            emiss_RWout				       <- mnl_RW_once(int1 = emiss_c_int[[i]][[q]][s,], Obs = cond_y[[s]][[q]][[i]], n_cat = q_emiss[q], mu_int_bar1 = c(t(emiss_mu_int_bar[[i]][[q]]) %*% xx[[1 + q]][s,]), V_int1 = emiss_V_int[[i]][[q]], scalar = emiss_scalar[[q]], candcov1 = emiss_candcov_comb)
            emiss[[s]][[q]][i,]		   <- PD_subj[[s]]$cat_emiss[iter, (sum(start[1:q]) + 1 + (i - 1) * q_emiss[q]):(sum(start[1:q]) + (i - 1) * q_emiss[q] + q_emiss[q])] <- (emiss_RWout$prob + .0001) / (1 + 0.0001 * q_emiss[q])
            emiss_naccept[[q]][s, i]	 <- emiss_naccept[[q]][s, i] + emiss_RWout$accept
            emiss_c_int[[i]][[q]][s,] <- emiss_RWout$draw_int
            emiss_int_subj[[s]][[q]][iter, (1 + (i - 1) * (q_emiss[q] - 1)) : ((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1))]	<- emiss_c_int[[i]][[q]][s,]
          }
        }
       } else if(data_distr == 'continuous'){
          # Sample subject values for normal emission distribution using Gibbs sampler   ---------

          # population level, conditional probabilities, separate for each dependent variable
          for(q in 1:n_dep){
            for (s in 1:n_subj){
              ss_subj[s] <- t(matrix(cond_y[[s]][[q]][[i]] - emiss_c_mu[[i]][[q]][s,1], nrow = 1) %*%
                                matrix(cond_y[[s]][[q]][[i]] - emiss_c_mu[[i]][[q]][s,1], ncol = 1))
              n_cond_y[s]       <- length(cond_y[[s]][[q]][[i]])
            }
            emiss_a_resvar_n <- sum(n_cond_y) / 2 + emiss_a0[[q]][i]
            emiss_b_resvar_n <- (sum(ss_subj) + 2 * emiss_b0[[q]][i]) / 2
            emiss_c_V[[i]][[q]] <- solve(stats::rgamma(1, shape = emiss_a_resvar_n, rate = emiss_b_resvar_n))
            emiss_sd_bar[[q]][iter, i] <- sqrt(emiss_c_V[[i]][[q]])
          }

          ### sampling subject specific means for the emission distributions, assuming known mean and var, see Lynch p. 244
          for(q in 1:n_dep){
            emiss_c_V_subj    <- (emiss_V_mu[[i]][[q]] * emiss_c_V[[i]][[q]]) / (2 * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
            for (s in 1:n_subj){
              emiss_mu0_subj_n  <- (emiss_V_mu[[i]][[q]] * sum(cond_y[[s]][[q]][[i]]) +  emiss_c_V[[i]][[q]] * c(t(emiss_c_mu_bar[[i]][[q]]) %*% xx[[q+1]][s,])) /
                (n_cond_y[s] * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
              emiss[[s]][[q]][i,1] <- PD_subj[[s]]$cont_emiss[iter, ((q - 1) * m + i)] <- emiss_c_mu[[i]][[q]][s,1] <- rnorm(1, emiss_mu0_subj_n, sqrt(emiss_c_V_subj))
              emiss[[s]][[q]][i,2] <- emiss_c_V[[i]][[q]]
              PD_subj[[s]]$cont_emiss[iter, (n_dep * m + (q - 1) * m + i)] <- sqrt(emiss_c_V[[i]][[q]])
            }
          }
       } else if(data_distr == 'count'){
         ### sampling subject specific means for the poisson emission distributions from a lognormal prior
         for(q in 1:n_dep){
           for (s in 1:n_subj){
             if(is.null(cond_y[[s]][[q]][[i]])){
               cond_y[[s]][[q]][[i]] <- round(exp(emiss_mu0_subj_bar),0)
             }
             # Sample subject specific lambda with RW-MH
             emiss_mu0_subj_bar <- c(t(emiss_c_mu_bar[[i]][[q]]) %*% xx[[q+1]][s,])
             emiss_candcov_comb <- (subj_data[[s]]$emiss_mhess[[q]][i] + emiss_V_mu[[i]][[q]]^-1)^-1
             emiss_rw_out <- pois_RW_once(lambda = emiss_c_mu[[i]][[q]][s,1],
                                          Obs = cond_y[[s]][[q]][[i]],
                                          mu_bar1 = emiss_mu0_subj_bar,
                                          V_1 = sqrt(emiss_V_mu[[i]][[q]]),
                                          scalar = emiss_scalar[[q]],
                                          candcov1 = emiss_candcov_comb)

             emiss[[s]][[q]][i,1]        <- PD_subj[[s]]$count_emiss[iter, ((q - 1) * m + i)] <- emiss_c_mu[[i]][[q]][s,1] <- emiss_rw_out$draw_lambda
             emiss_naccept[[q]][s, i]    <- emiss_naccept[[q]][s, i] + emiss_rw_out$accept
           }
         }
       }
    }



    # End of 1 MCMC iteration, save output values --------
    gamma_int_bar[iter, ]				   	<- unlist(lapply(gamma_mu_int_bar, "[",1,))
    if(nx[1] > 1){
      gamma_cov_bar[iter, ]      	<- unlist(lapply(gamma_mu_int_bar, "[",-1,))
    }
    gamma_V_int_bar[iter, ] <- unlist(lapply(gamma_V_int, function(e) as.vector(t(e))))
    gamma_prob_bar[iter,]			<- unlist(gamma_mu_prob_bar)
    if(data_distr == 'categorical'){
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
    } else if(data_distr == 'continuous'){
      for(q in 1:n_dep){
        emiss_mu_bar[[q]][iter, ]	<- as.vector(unlist(lapply(
          lapply(emiss_c_mu_bar, "[[", q), "[",1,)
        ))
        if(nx[1+q] > 1){
          emiss_cov_bar[[q]][iter, ]  <- as.vector(unlist(lapply(
            lapply(emiss_c_mu_bar, "[[", q), "[",-1,)
          ))
        }
        emiss_varmu_bar[[q]][iter,]	<- as.vector(unlist(sapply(emiss_V_mu, "[[", q)))
      }
    } else if(data_distr == 'count'){
      for(q in 1:n_dep){
        emiss_mu_bar[[q]][iter, ]	<- as.vector(unlist(lapply(
          lapply(emiss_c_mu_bar, "[[", q), "[",1,)
        ))
        if(nx[1+q] > 1){
          emiss_cov_bar[[q]][iter, ]  <- as.vector(unlist(lapply(
            lapply(emiss_c_mu_bar, "[[", q), "[",-1,)
          ))
        }
        emiss_varmu_bar[[q]][iter,]	<- as.vector(unlist(sapply(emiss_V_mu, "[[", q)))
      }
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
  if(data_distr == 'categorical'){
    if(return_path == TRUE){
      out <- list(input = list(data_distr = data_distr, m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,
                               burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj, gamma_int_subj = gamma_int_subj, emiss_int_subj = emiss_int_subj,
                  gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar, gamma_V_int_bar = gamma_V_int_bar,
                  emiss_int_bar = emiss_int_bar, emiss_cov_bar = emiss_cov_bar, gamma_prob_bar = gamma_prob_bar,
                  emiss_prob_bar = emiss_prob_bar, gamma_naccept = gamma_naccept, emiss_naccept = emiss_naccept,
                  sample_path = sample_path)
    } else {
      out <- list(input = list(data_distr = data_distr, m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,
                               burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj, gamma_int_subj = gamma_int_subj, emiss_int_subj = emiss_int_subj,
                  gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar, gamma_V_int_bar = gamma_V_int_bar,
                  emiss_int_bar = emiss_int_bar, emiss_cov_bar = emiss_cov_bar, gamma_prob_bar = gamma_prob_bar,
                  emiss_prob_bar = emiss_prob_bar, gamma_naccept = gamma_naccept, emiss_naccept = emiss_naccept)
    }
    class(out) <- append(class(out), c("mHMM", "cat"))
  } else if(data_distr == 'continuous'){
    if(return_path == TRUE){
      out <- list(input = list(data_distr = data_distr, m = m, n_dep = n_dep, J = J,
                               burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj, gamma_int_subj = gamma_int_subj,
                  gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar,
                  gamma_V_int_bar = gamma_V_int_bar,
                  emiss_cov_bar = emiss_cov_bar, gamma_prob_bar = gamma_prob_bar,
                  emiss_mu_bar = emiss_mu_bar, gamma_naccept = gamma_naccept,
                  emiss_varmu_bar = emiss_varmu_bar, emiss_sd_bar = emiss_sd_bar,
                  sample_path = sample_path)
    } else {
      out <- list(input = list(data_distr = data_distr, m = m, n_dep = n_dep, J = J,
                               burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj, gamma_int_subj = gamma_int_subj,
                  gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar,
                  gamma_V_int_bar = gamma_V_int_bar,
                  emiss_cov_bar = emiss_cov_bar, gamma_prob_bar = gamma_prob_bar,
                  emiss_mu_bar = emiss_mu_bar, gamma_naccept = gamma_naccept,
                  emiss_varmu_bar = emiss_varmu_bar, emiss_sd_bar = emiss_sd_bar)
    }
    class(out) <- append(class(out), c("mHMM", "cont"))
  } else if(data_distr == 'count'){
    if(return_path == TRUE){
      out <- list(input = list(data_distr = data_distr, m = m, n_dep = n_dep, J = J,
                               burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj,
                  gamma_int_subj = gamma_int_subj,
                  gamma_int_bar = gamma_int_bar,
                  gamma_cov_bar = gamma_cov_bar,
                  gamma_V_int_bar = gamma_V_int_bar,
                  emiss_cov_bar = emiss_cov_bar,
                  gamma_prob_bar = gamma_prob_bar,
                  emiss_mu_bar = emiss_mu_bar,
                  emiss_varmu_bar = emiss_varmu_bar,
                  gamma_naccept = gamma_naccept,
                  emiss_naccept = emiss_naccept,
                  sample_path = sample_path)
    } else {
      out <- list(input = list(data_distr = data_distr, m = m, n_dep = n_dep, J = J,
                               burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj,
                  gamma_int_subj = gamma_int_subj,
                  gamma_int_bar = gamma_int_bar,
                  gamma_cov_bar = gamma_cov_bar,
                  gamma_V_int_bar = gamma_V_int_bar,
                  gamma_prob_bar = gamma_prob_bar,
                  emiss_mu_bar = emiss_mu_bar,
                  emiss_cov_bar = emiss_cov_bar,
                  emiss_varmu_bar = emiss_varmu_bar,
                  gamma_naccept = gamma_naccept,
                  emiss_naccept = emiss_naccept)
    }
    class(out) <- append(class(out), c("mHMM","count"))
  }
  return(out)
}
