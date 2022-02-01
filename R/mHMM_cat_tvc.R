#' Multilevel hidden  Markov model with time varying covariate to predict gamma using Bayesian estimation
#'
#' \code{mHMM_cat_tv} fits a multilevel (also known as mixed or random effects)
#' hidden Markov model (HMM) with ONE time varying covariate to predict gamma.
#' The model is designed to be used on intense longitudinal data with
#' categorical observations of multiple subjects and is fitted using Bayesian
#' estimation. The function \code{mHMM_cat_tv} creates an object of class mHMM.
#' By using a multilevel framework, we allow for heterogeneity in the model
#' parameters between subjects, while estimating one overall HMM. In addition to
#' one time varying covariate, the function includes the possibility to add time
#' invariant covariates at level 2 (i.e., at the subject level) and have varying
#' observation lengths over subjects. For a short description of the package see
#' \link{mHMMbayes}. See \code{vignette("tutorial-mhmm")} for an introduction to
#' multilevel hidden Markov models and the package, and see
#' \code{vignette("estimation-mhmm")} for an overview of the used estimation
#' algorithms.
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
#' If covariates are specified and the user wants to set the values for the
#' parameters of the hyper-prior distributions manually, the specification of
#' the elements in the arguments of the hyper-prior parameter values of the
#' normal distribution on the means change as follows: the number of rows in the
#' matrices \code{gamma_mu0_int} and \code{emiss_mu0} are equal to 1 + the number of
#' covariates used to predict the transition probability matrix for
#' \code{gamma_mu0_int} and the emission distribution for \code{emiss_mu0} (i.e.,
#' the first row correspond to the hyper-prior mean values of the intercepts,
#' the subsequent rows correspond to the hyper-prior mean values of the
#' regression coefficients connected to each of the covariates), and
#' \code{gamma_K0_int} and \code{emiss_K0} are now a matrix with the number of
#' hypothetical prior subjects on which the vectors of means (i.e., the rows in
#' \code{gamma_mu0_int} or \code{emiss_mu0}) are based on the diagonal, and
#' off-diagonal elements equal to 0. Note that the hyper-prior parameter values
#' of the inverse Wishart distribution on the covariance matrix remains
#' unchanged, as the estimates of the regression coefficients for the covariates
#' are fixed over subjects.
#'
#' @param s_data A matrix containing the observations to be modelled, where the
#'   rows represent the observations over time. In \code{s_data}, the first
#'   column indicates subject id number. Hence, the id number is repeated over
#'   rows equal to the number of observations for that subject. The subsequent
#'   \code{n_dep} columns contain the dependent variable(s). The last column
#'   contains the time varying covariate used to predict the transition
#'   probability matrix gamma. Note that the dependent variables have to be
#'   numeric, i.e., they cannot be a (set of) factor variable(s). The total
#'   number of rows are equal to the sum over the number of observations of each
#'   subject, and the number of columns are equal to the number of dependent
#'   variables (\code{n_dep}) + 2. The number of observations can vary over
#'   subjects.
#' @param gen List containing the following elements denoting the general model
#'   properties:
#'   \itemize{\item{\code{m}: numeric vector with length 1 denoting the number
#'   of hidden states}
#'   \item{\code{n_dep}: numeric vector with length 1 denoting the
#'   number of dependent variables}
#'   \item{\code{q_emiss}: numeric vector with length \code{n_dep} denoting the
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
#'   values for gamma. The subsequent elements contain \code{m} by
#'   \code{q_emiss[k]} matrices for the start values for each of the \code{k} in
#'   \code{n_dep} emission distribution(s). Note that \code{start_val} should not contain
#'   nested lists (i.e., lists within lists).
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
#' @param print_iter The argument print_iter is depricated; please use
#'   show_progress instead to show the progress of the algorithm.
#' @param show_progress A logical scaler. Should the function show a text
#'   progress bar in the \code{R} console to represent the progress of the
#'   algorithm (\code{show_progress = TRUE}) or not (\code{show_progress =
#'   FALSE}). Defaults to \code{show_progress = TRUE}.
#' @param gamma_hyp_prior An optional list containing user specified parameters
#'  of the hyper-prior distribution on the multivariate normal distribution
#'  of the intercepts (and regression coefficients given that covariates are
#'  used) of the multinomial regression model of the transition probability
#'  matrix gamma. The hyper-prior of the mean intercepts is a multivariate
#'  Normal distribution, the hyper-prior of the covariance matrix between the
#'  set of (state specific) intercepts is an Inverse Wishart distribution.
#'
#'  Hence, the list \code{gamma_hyp_prior} contains the following elements:
#'  \itemize{\item{\code{gamma_mu0_int}: a list containing m matrices; one matrix
#'  for each row of the transition probability matrix gamma. Each matrix
#'  contains the hypothesized mean values of the intercepts. Hence, each matrix
#'  consists of one row (when not including covariates in the model) and
#'  \code{m} - 1 columns}
#'  \item{\code{gamma_K0_int}: numeric vector with length 1 denoting the number of
#'  hypothetical prior subjects on which the vector of means \code{gamma_mu0_int} is
#'  based}
#'  \item{\code{gamma_nu_int}: numeric vector with length 1 denoting the degrees of
#'  freedom of the Inverse Wishart distribution}
#'  \item{\code{gamma_V0_int}: matrix of \code{m} - 1 by \code{m} - 1 containing the
#'  hypothesized variance-covariance matrix between the set of intercepts.}}
#'  Note that \code{gamma_K0_int}, \code{gamma_nu_int} and \code{gamma_V0_int} are assumed
#'  equal over the states. The mean values of the intercepts (and regression
#'  coefficients of the covariates) denoted by \code{gamma_mu0_int} are allowed to
#'  vary over the states.
#'
#'  The default values for the hyper-prior on gamma are: all elements of the
#'  matrices contained in \code{gamma_mu0_int} set to 0, \code{gamma_K0_int} set to 1,
#'  \code{gamma_nu_int} set to 3 + m - 1, and the diagonal of \code{gamma_V0_int} (i.e.,
#'  the variance) set to 3 + m - 1 and the off-diagonal elements (i.e., the
#'  covariance) set to 0.
#'
#'  See \emph{Details} below if covariates are used for changes in the settings
#'  of the arguments of \code{gamma_hyp_prior}.
#' @param emiss_hyp_prior An optional list containing user specified parameters
#'  of the hyper-prior distribution on the multivariate normal distribution
#'  of the intercepts (and regression coefficients given that covariates are
#'  used) of the multinomial regression model of the emission distribution.
#'  The hyper-prior for the mean intercepts is a multivariate
#'  Normal distribution, the hyper-prior for the covariance matrix between the
#'  set of (state specific) intercepts is an Inverse Wishart distribution.
#'
#'  Hence, the list \code{emiss_hyp_prior} contains the following elements:
#'  \itemize{\item{\code{emiss_mu0}: a list of lists: \code{emiss_mu0} contains
#'  \code{n_dep} lists, i.e., one list for each dependent variable \code{k}. Each
#'  of these lists contains m matrices; one matrix for each set of emission
#'  probabilities within a state. The matrices contain the hypothesized mean
#'  values of the intercepts. Hence, each matrix consists of one row (when not
#'  including covariates in the model) and \code{q_emiss[k]} - 1 columns}
#'  \item{\code{emiss_K0}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables, where each element is a numeric vector
#'  with length 1 denoting the number of hypothetical prior subjects on which
#'  the vector of means \code{emiss_mu0} is based}
#'  \item{\code{emiss_nu}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables, where each element is a numeric vector
#'  with length 1 denoting the degrees of freedom of the Inverse Wishart
#'  distribution}
#'  \item{\code{emiss_V}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables \code{k}, where each element is a matrix
#'  of \code{q_emiss[k]} - 1 by \code{q_emiss[k]} - 1 containing the
#'  hypothesized variance-covariance matrix between the set of intercepts.}}
#'  Note that \code{emiss_K0}, \code{emiss_nu} and \code{emiss_V} are assumed
#'  equal over the states. The mean values of the intercepts (and regression
#'  coefficients of the covariates) denoted by \code{emiss_mu0} are allowed to
#'  vary over the states.
#'
#'  The default values for the hyper-prior on the emission distribution(s) are:
#'  all elements of the matrices contained in \code{emiss_mu0} set to 0,
#'  \code{emiss_K0} set to 1, \code{emiss_nu} set to 3 + \code{q_emiss[k]} - 1,
#'  and the diagonal of \code{gamma_V0_int} (i.e., the variance) set to 3 +
#'  \code{q_emiss[k]} - 1 and the off-diagonal elements (i.e., the covariance)
#'  set to 0.
#'
#'  See \emph{Details} below if covariates are used for changes in the settings
#'  of the arguments of \code{emiss_hyp_prior}.
#' @param gamma_sampler An optional list containing user specified settings for
#'   the proposal distribution of the random walk (RW) Metropolis sampler for
#'   the subject level parameter estimates of the intercepts modeling the
#'   transition probability matrix. The list \code{gamma_sampler} contains the
#'   following elements:
#'  \itemize{\item{\code{gamma_int_mle0}: a numeric vector with length \code{m}
#'  - 1 denoting the start values for the maximum likelihood estimates of the
#'  intercepts for the transition probability matrix gamma, based on the pooled
#'  data (data over all subjects)}
#'  \item{\code{gamma_scalar}: a numeric vector with length 1 denoting the scale
#'  factor \code{s}. That is, The scale of the proposal distribution is composed
#'  of a covariance matrix Sigma, which is then tuned by multiplying it by a
#'  scaling factor \code{s}^2}
#'  \item{\code{gamma_w}: a numeric vector with length 1 denoting the weight for
#'  the overall log likelihood (i.e., log likelihood based on the pooled data
#'  over all subjects) in the fractional likelihood.}}
#'   Default settings are: all elements in \code{gamma_int_mle0} set to 0,
#'   \code{gamma_scalar} set to 2.93 / sqrt(\code{m} - 1), and \code{gamma_w} set to
#'   0.1. See the section \emph{Scaling the proposal distribution of the RW
#'   Metropolis sampler} in \code{vignette("estimation-mhmm")} for details.
#' @param emiss_sampler An optional list containing user specified settings for
#'   the proposal distribution of the random walk (RW) Metropolis sampler for
#'   the subject level parameter estimates of the intercepts modeling the
#'   emission distributions of the dependent variables \code{k}. The list
#'   \code{emiss_sampler} contains the following elements:
#'  \itemize{\item{\code{emiss_int_mle0}: a list containing \code{n_dep} elements
#'  corresponding to each of the dependent variables \code{k}, where each
#'  element is a a numeric vector with length \code{q_emiss[k]} - 1 denoting the
#'  start values for the maximum likelihood estimates of the intercepts for
#'  the emission distribution, based on the pooled data (data over all
#'  subjects)}
#'  \item{\code{emiss_scalar}: a list containing \code{n_dep} elements
#'  corresponding to each of the dependent variables, where each element is a
#'  numeric vector with length 1 denoting the scale factor \code{s}. That is,
#'  The scale of the proposal distribution is composed of a covariance matrix
#'  Sigma, which is then tuned by multiplying it by a scaling factor \code{s}^2}
#'  \item{\code{emiss_w}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables, where each element is a numeric vector
#'  with length 1 denoting the weight for the overall log likelihood (i.e., log
#'  likelihood based on the pooled data over all subjects) in the fractional
#'  likelihood.}}
#'   Default settings are: all elements in \code{emiss_int_mle0} set to 0,
#'   \code{emiss_scalar} set to 2.93 / sqrt(\code{q_emiss[k]} - 1), and
#'   \code{emiss_w} set to 0.1. See the section \emph{Scaling the proposal
#'   distribution of the RW Metropolis sampler} in
#'   \code{vignette("estimation-mhmm")} for details.
#'
#' @return \code{mHMM} returns an object of class \code{mHMM}, which has
#'   \code{print} and \code{summary} methods to see the results.
#'   The object contains the following components:
#'   \describe{
#'   \item{\code{PD_subj}}{A list containing one matrix per subject with the
#'   subject level parameter estimates and the log likelihood over the
#'   iterations of the hybrid Metropolis within Gibbs sampler. The iterations of
#'   the sampler are contained in the rows, and the columns contain the subject
#'   level (parameter) estimates of subsequently the emission probabilities, the
#'   transition probabilities and the log likelihood.}
#'   \item{\code{gamma_prob_bar}}{A matrix containing the group level parameter
#'   estimates of the transition probabilities over the iterations of the hybrid
#'   Metropolis within Gibbs sampler. The iterations of the sampler are
#'   contained in the rows, and the columns contain the group level parameter
#'   estimates. If covariates were included in the analysis, the group level
#'   probabilities represent the predicted probability given that the covariate
#'   is at the average value for continuous covariates, or given that the
#'   covariate equals zero for dichotomous covariates.}
#'   \item{\code{gamma_int_bar}}{A matrix containing the group level intercepts
#'   of the multinomial logistic regression modeling the transition
#'   probabilities over the iterations of the hybrid Metropolis within Gibbs
#'   sampler. The iterations of the sampler are contained in the rows, and the
#'   columns contain the group level intercepts.}
#'   \item{\code{gamma_cov_bar}}{A matrix containing the group level regression
#'   coefficients of the multinomial logistic regression predicting the
#'   transition probabilities over the iterations of the hybrid Metropolis within
#'   Gibbs sampler. The iterations of the sampler are contained in the rows, and
#'   the columns contain the group level regression coefficients.}
#'   \item{\code{gamma_int_subj}}{A list containing one matrix per subject
#'   denoting the subject level intercepts of the multinomial logistic
#'   regression modeling the transition probabilities over the iterations of the
#'   hybrid Metropolis within Gibbs sampler. The iterations of the sampler are
#'   contained in the rows, and the columns contain the subject level
#'   intercepts.}
#'   \item{\code{gamma_int_naccept}}{A matrix containing the number of accepted
#'   draws at the subject level RW Metropolis step for each set of parameters of
#'   the transition probabilities. The subjects are contained in the rows, and
#'   the columns contain the sets of parameters.}
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
#'   the multinomial logistic regression modeling the probabilities of the
#'   emission distribution over the iterations of the hybrid Metropolis within
#'   Gibbs sampler. The iterations of the sampler are contained in the rows of
#'   the matrix, and the columns contain the group level intercepts.}
#'   \item{\code{emiss_cov_bar}}{A list containing one matrix per dependent
#'   variable, denoting the group level regression coefficients of the
#'   multinomial logistic regression predicting the emission probabilities within
#'   each of the dependent variables over the iterations of the hybrid
#'   Metropolis within Gibbs sampler. The iterations of the sampler are
#'   contained in the rows  of the matrix, and the columns contain the group
#'   level regression coefficients.}
#'   \item{\code{emiss_int_subj}}{A list containing one list per subject denoting
#'   the subject level intercepts of each dependent variable of the multinomial
#'   logistic regression modeling the probabilities of the emission distribution
#'   over the iterations of the hybrid Metropolis within Gibbs sampler. Each
#'   lower level list contains one matrix per dependent variable, in which
#'   iterations of the sampler are contained in the rows, and the columns
#'   contain the subject level intercepts.}
#'   \item{\code{emiss_naccept}}{A list containing one matrix per dependent
#'   variable with the number of accepted draws at the subject level RW
#'   Metropolis step for each set of parameters of the emission distribution.
#'   The subjects are contained in the rows, and the columns of the matrix
#'   contain the sets of parameters.}
#'   \item{\code{input}}{Overview of used input specifications: the number of
#'   states \code{m}, the number of used dependent variables \code{n_dep}, the
#'   number of output categories for each of the dependent variables
#'   \code{q_emiss}, the number of iterations \code{J} and the specified burn in
#'   period \code{burn_in} of the hybrid Metropolis within Gibbs sampler, the
#'   number of subjects \code{n_subj}, the observation length for each subject
#'   \code{n_vary}, and the column names of the dependent variables
#'   \code{dep_labels}.}
#'   \item{\code{sample_path}}{A list containing one matrix per subject with the
#'   sampled hidden state sequence over the hybrid Metropolis within Gibbs
#'   sampler. The time points of the dataset are contained in the rows, and the
#'   sampled paths over the iterations are contained in the columns. Only
#'   returned if \code{return_path = TRUE}. }
#' }
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
#' ###### Example on package example data
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
#' # Run a model including a covariate. Here, the covariate (standardized CDI
#' # change) predicts the emission distribution for each of the 4 dependent
#' # variables:
#'
#' n_subj <- 10
#' xx <- rep(list(matrix(1, ncol = 1, nrow = n_subj)), (n_dep + 1))
#' for(i in 2:(n_dep + 1)){
#'   xx[[i]] <- cbind(xx[[i]], nonverbal_cov$std_CDI_change)
#' }
#' out_2st_c <- mHMM(s_data = nonverbal, xx = xx,
#'                   gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                   start_val = c(list(start_TM), start_EM),
#'                   mcmc = list(J = 11, burn_in = 5))
#'
#' }
#' ###### Example on simulated data
#' # Simulate data for 10 subjects with each 100 observations:
#' n_t <- 100
#' n <- 10
#' m <- 2
#' q_emiss <- 3
#' gamma <- matrix(c(0.8, 0.2,
#'                   0.3, 0.7), ncol = m, byrow = TRUE)
#' emiss_distr <- list(matrix(c(0.5, 0.5, 0.0,
#'                         0.1, 0.1, 0.8), nrow = m, ncol = q_emiss, byrow = TRUE))
#' data1 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
#'                   emiss_distr = emiss_distr, var_gamma = .5, var_emiss = .5)
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
#'
#' @export
#'
#'

mHMM_cat_tv <- function(s_data, gen, xx = NULL, start_val, mcmc, return_path = FALSE, print_iter, show_progress = TRUE,
                 gamma_hyp_prior = NULL, emiss_hyp_prior = NULL, gamma_sampler = NULL, emiss_sampler = NULL){

  if(!missing(print_iter)){
    warning("The argument print_iter is depricated; please use show_progress instead to show the progress of the algorithm.")
  }
  # Initialize data -----------------------------------
  # dependent variable(s), sample size, dimensions gamma and conditional distribuiton
  n_dep			 <- gen$n_dep
  if(ncol(s_data) != n_dep + 2){
    stop(paste("The number of columns in s_data should equal ", n_dep + 2, ":", "\n", "first column denoting subject id, subsequently ", n_dep, " column(s) for the dependent variable(s), and a final column denoting the time varying covariate.", sep = ""))
  }
  dep_labels <- colnames(s_data[,2:(n_dep+1)])
  id         <- unique(s_data[,1])
  n_subj     <- length(id)
  subj_data  <- rep(list(NULL), n_subj)
  if(sum(sapply(s_data, is.factor)) > 0 ){
    stop("Your data contains factorial variables, which cannot be used as input in the function mHMM. All variables have to be numerical.")
  }
  for(s in 1:n_subj){
    subj_data[[s]]$y <- as.matrix(s_data[s_data[,1] == id[s],][,2:(n_dep + 1)], ncol = n_dep)
    subj_data[[s]]$xx_t <- matrix(c(rep(1, sum(s_data[,1] == id[s])), s_data[s_data[,1] == id[s],][, n_dep + 2]), ncol = 2)
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
    subj_data[[s]]	<- c(subj_data[[s]], n = n, list(gamma_converge = numeric(m), gamma_int_mle = matrix(, m, (m - 1) * 2),
                                                    gamma_mhess = matrix(, (m - 1) * 2, (m - 1) * 2), emiss_converge =
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
    if(!is.list(xx) | length(xx) != n_dep1){
      stop("If xx is specified, xx should be a list, with the number of elements equal to the number of dependent variables + 1")
    }
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
    gamma_int_mle0  <- rep(0, m - 1)
    gamma_bet_mle0  <- rep(0, m - 1)
    gamma_scalar    <- 2.93 / sqrt(m - 1)
    gamma_w         <- .1
  } else {
    gamma_int_mle0  <- gamma_sampler$gamma_int_mle0
    gamma_bet_mle0  <- gamma_sampler$gamma_bet_mle0
    gamma_scalar    <- gamma_sampler$gamma_scalar
    gamma_w         <- gamma_sampler$gamma_w
  }

  # Initialize emiss sampler
  if(is.null(emiss_sampler)){
    emiss_int_mle0 <- rep(list(NULL), n_dep)
    emiss_scalar	<- rep(list(NULL), n_dep)
    for(q in 1:n_dep){
      emiss_int_mle0[[q]] <- rep(0, q_emiss[q] - 1)
      emiss_scalar[[q]] 	<- 2.93 / sqrt(q_emiss[q] - 1)
    }
    emiss_w		<- .1
  } else {
    emiss_int_mle0	<- emiss_sampler$emiss_int_mle0
    emiss_scalar 	<- emiss_sampler$emiss_scalar
    emiss_w    		<- emiss_sampler$emiss_w
  }

  # Initialize Gamma hyper prior
  if(is.null(gamma_hyp_prior)){
    gamma_mu0_int	    <- rep(list(matrix(0,nrow = nx[1], ncol = m - 1)), m)
    gamma_K0_int			<- diag(1, nx[1])
    gamma_nu_int			<- 3 + m - 1
    gamma_V0_int			<- gamma_nu_int * diag(m - 1)
    gamma_mu0_bet     <- rep(0, m-1)
    gamma_V0_bet      <- (3 + m - 1) * diag(m-1)
  } else {
    ###### BUILD in a warning / check if gamma_mu0_int is a matrix when given, with  nrows equal to the number of covariates
    gamma_mu0_int			<- gamma_hyp_prior$gamma_mu0_int
    gamma_K0_int			<- gamma_hyp_prior$gamma_K0_int
    gamma_nu_int			<- gamma_hyp_prior$gamma_nu_int
    gamma_V0_int			<- gamma_hyp_prior$gamma_V0_int
    gamma_mu0_bet     <- gamma_hyp_prior$gamma_mu0_bet
    gamma_V0_bet      <- gamma_hyp_prior$gamma_V0_bet
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
      emiss_mu0[[q]]	 <- emiss_hyp_prior$emiss_mu0[[q]]
      emiss_nu[[q]]	 <- emiss_hyp_prior$emiss_nu[[q]]
      emiss_V[[q]]		 <- emiss_hyp_prior$emiss_V[[q]]
      emiss_K0[[q]]	 <- diag(emiss_hyp_prior$emiss_K0, nx[1 + q])
    }
  }


  # Define objects used to store data in mcmc algorithm, not returned ----------------------------
  # overall
  c <- llk <- numeric(1)
  sample_path <- lapply(n_vary, dif_matrix, cols = J)
  trans <- xx_t_state <- rep(list(vector("list", m)), n_subj)


  # gamma
  gamma_int_mle_pooled <- gamma_pooled_ll <- gamma_hess_mle_pooled <- vector("list", m)
  gamma_c_int <- rep(list(matrix(, n_subj, (m-1))), m)
  gamma_c_bet <- rep(list(rep(0, m-1)),m)
  gamma_mu_int_bar <- gamma_V_int <- vector("list", m)
  gamma_mu_prob_bar <- rep(list(numeric(m)), m)
  gamma_int_naccept <- matrix(0, n_subj, m)
  gamma_bet_naccept <- matrix(0, n_subj, m)

  # emiss
  cond_y <- lapply(rep(n_dep, n_subj), nested_list, m = m)
  emiss_int_mle_pooled <- emiss_pooled_ll <- rep(list(vector("list", n_dep)), m)
  emiss_c_int <- rep(list(lapply(q_emiss - 1, dif_matrix, rows = n_subj)), m)
  emiss_mu_int_bar <- emiss_V_int <- rep(list(vector("list", n_dep)), m)
  emiss_mu_prob_bar <- rep(list(lapply(q_emiss, dif_vector)), m)
  emiss_naccept <- rep(list(matrix(0, n_subj, m)), n_dep)


  # Define objects that are returned from mcmc algorithm ----------------------------
  # Define object for subject specific posterior density, put start values on first row
  if(length(start_val) != n_dep + 1){
    stop("The number of elements in the list start_val should be equal to 1 + the number of dependent variables,
         and should not contain nested lists (i.e., lists within lists)")
  }
  PD 					  <- matrix(, nrow = J, ncol = sum(m * q_emiss) + m * m + 1)
  PD_emiss_names   <- paste("q", 1, "_emiss", rep(1:q_emiss[1], m), "_S", rep(1:m, each = q_emiss[1]), sep = "")
  if(n_dep > 1){
    for(q in 2:n_dep){
      PD_emiss_names <- c(PD_emiss_names, paste("q", q, "_emiss", rep(1:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q]), sep = ""))
    }
  }
  colnames(PD) 	<- c(PD_emiss_names, paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = ""), "LL")

  PD[1, ((sum(m * q_emiss) + 1)) :((sum(m * q_emiss) + m * m))] <- unlist(sapply(start_val, t))[1:(m*m)]
  PD[1, 1:((sum(m * q_emiss)))] <- unlist(lapply(start_val, t))[(m*m + 1): (m*m + sum(m * q_emiss))]

  PD_subj				<- rep(list(PD), n_subj)

  # Define object for population posterior density (probabilities and regression coefficients parameterization )
  gamma_prob_bar		<- matrix(, nrow = J, ncol = (m * m))
  colnames(gamma_prob_bar) <- paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = "")
  gamma_prob_bar[1,] <- PD[1,(sum(m*q_emiss) + 1):(sum(m * q_emiss) + m * m)]
  emiss_prob_bar			<- lapply(q_emiss * m, dif_matrix, rows = J)
  names(emiss_prob_bar) <- dep_labels
  for(q in 1:n_dep){
    colnames(emiss_prob_bar[[q]]) <- paste("Emiss", rep(1:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q]), sep = "")
    start <- c(0, q_emiss * m)
    emiss_prob_bar[[q]][1,] <- PD[1,(sum(start[1:q]) + 1):(sum(start[1:q]) + (m * q_emiss[q]))]
  }
  gamma_int_bar				<- matrix(, nrow = J, ncol = ((m-1) * m))
  colnames(gamma_int_bar) <- paste("int_S", rep(1:m, each = m-1), "toS", rep(2:m, m), sep = "")
  gamma_int_bar[1,] <- as.vector(prob_to_int(matrix(gamma_prob_bar[1,], byrow = TRUE, ncol = m, nrow = m)))
  if(nx[1] > 1){
    gamma_cov_bar				<- matrix(, nrow = J, ncol = ((m-1) * m) * (nx[1] - 1))
    colnames(gamma_cov_bar) <- paste( paste("cov", 1 : (nx[1] - 1), "_", sep = ""), "S", rep(1:m, each = (m-1) * (nx[1] - 1)), "toS", rep(2:m, m * (nx[1] - 1)), sep = "")
    gamma_cov_bar[1,] <- 0
  } else{
    gamma_cov_bar <- "No covariates where used to predict the transition probability matrix"
  }
  gamma_tv_cov         <- matrix(, ncol = (m-1)*m, nrow = J)
  colnames(gamma_tv_cov) <- paste("tv_cov","_S", rep(1:m, each = (m-1)), "toS", rep(2:m, m), sep = "")
  emiss_int_bar			<- lapply((q_emiss-1) * m, dif_matrix, rows = J)
  names(emiss_int_bar) <- dep_labels
  for(q in 1:n_dep){
    colnames(emiss_int_bar[[q]]) <-  paste("int_Emiss", rep(2:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q] - 1), sep = "")
    emiss_int_bar[[q]][1,] <- as.vector(prob_to_int(matrix(emiss_prob_bar[[q]][1,], byrow = TRUE, ncol = q_emiss[q], nrow = m)))
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
  if(show_progress == TRUE){
    cat("Progress of the Bayesian mHMM algorithm:", "\n")
    pb <- utils::txtProgressBar(min = 2, max = J, style = 3)
  }
  for (iter in 2 : J){

    # For each subject, obtain sampled state sequence with subject individual parameters ----------
    for(s in 1:n_subj){
      # Run forward algorithm, obtain subject specific forward proababilities and log likelihood
      gamma_t       <- timevar_gamma(int = matrix(gamma_int_subj[[s]][iter-1,], byrow = TRUE, ncol = m-1),
                                     bet = matrix(unlist(gamma_c_bet), byrow = TRUE, ncol = m-1),
                                     xx_t = subj_data[[s]]$xx_t, m = m)
      forward				<- cat_tv_mult_fw_r_to_cpp(x = subj_data[[s]]$y, m = m, emiss = emiss[[s]], tgamma = gamma_t, n_dep = n_dep, delta=NULL)
      alpha         <- forward[[1]]
      c             <- max(forward[[2]][, subj_data[[s]]$n])
      llk           <- c + log(sum(exp(forward[[2]][, subj_data[[s]]$n] - c)))
      PD_subj[[s]][iter, sum(m * q_emiss) + m * m + 1] <- llk

      # Using the forward probabilities, sample the state sequence in a backward manner.
      # In addition, saves:
        # state transitions in trans,
        # conditional observations within states in cond_y
        # observations of the time dependent covariate per sampled state in xx_t_state
      trans[[s]]					                  <- vector("list", m)
      xx_t_state[[s]]                       <- vector("list", m)
      sample_path[[s]][n_vary[[s]], iter] 	<- sample(1:m, 1, prob = c(alpha[, n_vary[[s]]]))
      for(t in (subj_data[[s]]$n - 1):1){
        sample_path[[s]][t,iter] 	              <- sample(1:m, 1, prob = (alpha[, t] * matrix(gamma_t[t,], byrow = TRUE, ncol = m)[,sample_path[[s]][t + 1, iter]]))
        trans[[s]][[sample_path[[s]][t,iter]]]	<- c(trans[[s]][[sample_path[[s]][t, iter]]], sample_path[[s]][t + 1, iter])
      }
      for (i in 1:m){
        trans[[s]][[i]] <- rev(trans[[s]][[i]])
        xx_t_state[[s]][[i]] <- subj_data[[s]]$xx_t[-n_vary[s],][which(sample_path[[s]][-n_vary[s],iter] == i),2]
        for(q in 1:n_dep){
          cond_y[[s]][[i]][[q]] <- subj_data[[s]]$y[sample_path[[s]][, iter] == i, q]
        }
      }
    }

    # The remainder of the mcmc algorithm is state specific
    for(i in 1:m){

      # Obtain pooled over all subjects for optim and RW for beta
      Obs_p               <- matrix(c(rep(1:n_subj, times = sapply(trans, lengths)[i,]),
                                      unlist(sapply(trans, "[[", i))), ncol = 2)
      xx_t_p              <- matrix(c(rep(1, sum(sapply(xx_t_state, lengths)[i,])),
                                      unlist(sapply(xx_t_state, "[[", i))), ncol = 2)

      # Obtain MLE of the covariance matrices and log likelihood of gamma and emiss at subject and population level -----------------
      # used to scale the propasal distribution of the RW Metropolis sampler

      # population level, transition matrix
      X                   <- makeX(xx_t = matrix(c(rep(1, length(xx_t_p[,2]) + m), c(xx_t_p[,2], rep(0,m))), ncol = 2), n_cat = m)
      gamma_mle_pooled		<- optim(c(gamma_int_mle0, gamma_bet_mle0), llmnl_bet_optim, Obs = c(Obs_p[,2], 1:m),
                                 X = X, n_cat = m, method = "BFGS", hessian = FALSE,
                                 control = list(fnscale = -1))
      gamma_int_mle_pooled[[i]]  <- gamma_mle_pooled$par
      gamma_hess_mle_pooled[[i]]   <- mnlHess_bet(bet = gamma_mle_pooled$par, n_cat = m,
                                             Obs = c(Obs_p[,2], 1:m),
                                             xx_t = rbind(xx_t_p, matrix(rep(c(0,1), each = m), ncol = 2)))
      gamma_pooled_ll[[i]]			<- gamma_mle_pooled$value

      # population level, conditional probabilities, seperate for each dependent variable
      for(q in 1:n_dep){
        cond_y_pooled					      <- numeric()
        ### MOET OOK ECHT BETER KUNNEN, eerst # cond_y_pooled				<- unlist(sapply(cond_y, "[[", m))
        for(s in 1:n_subj){
          cond_y_pooled             <- c(cond_y_pooled, cond_y[[s]][[i]][[q]])
        }
        emiss_mle_pooled		<- optim(emiss_int_mle0[[q]], llmnl_int, Obs = c(cond_y_pooled, c(1:q_emiss[q])),
                                   n_cat = q_emiss[q], method = "BFGS", hessian = TRUE,
                                   control = list(fnscale = -1))
        emiss_int_mle_pooled[[i]][[q]]  <- emiss_mle_pooled$par
        emiss_pooled_ll[[i]][[q]]				<- emiss_mle_pooled$value
      }

      # subject level
      for (s in 1:n_subj){
        wgt 				<- subj_data[[s]]$n / n_total

        # subject level, transition matrix
        X                 <- makeX(xx_t = matrix(c(rep(1, length(xx_t_state[[s]][[i]]) + m), c(xx_t_state[[s]][[i]], rep(0,m))), ncol = 2), n_cat = m)
        gamma_out					<- optim(gamma_int_mle_pooled[[i]], llmnl_bet_frac, Obs = c(trans[[s]][[i]], c(1:m)),
                               X = X, n_cat = m, pooled_likel = gamma_pooled_ll[[i]], w = gamma_w, wgt = wgt,
                               method="BFGS", hessian = FALSE, control = list(fnscale = -1))
        if(gamma_out$convergence == 0){
          subj_data[[s]]$gamma_converge[i] <- 1
          subj_data[[s]]$gamma_int_mle[i,] <- gamma_out$par
          subj_data[[s]]$gamma_mhess	<- mnlHess_bet(bet = gamma_out$par, n_cat = m,
                                                    Obs = c(trans[[s]][[i]], c(1:m)),
                                                    xx_t = matrix(c(rep(1, length(xx_t_state[[s]][[i]]) + m),
                                                                    c(xx_t_state[[s]][[i]], rep(0,m))), ncol = 2))
        } else {
          subj_data[[s]]$gamma_converge[i] <- 0
          subj_data[[s]]$gamma_int_mle[i,] <- rep(0, m - 1)
          subj_data[[s]]$gamma_mhess	<- diag((m-1)*2)
        }
        # if this is first iteration, use MLE for current values RW metropolis sampler
        if (iter == 2){
          gamma_c_int[[i]][s,]		<- gamma_out$par[1:(m-1)]
        }

        # subject level, conditional probabilities, seperate for each dependent variable
        for(q in 1:n_dep){
          emiss_out				<- optim(emiss_int_mle_pooled[[i]][[q]], llmnl_int_frac, Obs = c(cond_y[[s]][[i]][[q]], c(1:q_emiss[q])),
                                n_cat = q_emiss[q], pooled_likel = emiss_pooled_ll[[i]][[q]],
                                w = emiss_w, wgt = wgt, method = "BFGS", hessian = TRUE, control = list(fnscale = -1))
          if(emiss_out$convergence == 0){
            subj_data[[s]]$emiss_converge[[q]][i]	 <- 1
            subj_data[[s]]$emiss_int_mle[[q]][i,] <- emiss_out$par
            subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ]		<-
              mnlHess_int(int = subj_data[[s]]$emiss_int_mle[[q]][i,], Obs = c(cond_y[[s]][[i]][[q]], c(1:q_emiss[q])), n_cat =  q_emiss[q])
          } else {
            subj_data[[s]]$emiss_converge[[q]][i]	 <- 0
            subj_data[[s]]$emiss_int_mle[[q]][i,]  <- rep(0, q_emiss[q] - 1)
            subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ]	<- diag(q_emiss[q] - 1)
          }
          # if this is first iteration, use MLE for current values RW metropolis sampler
          if (iter == 2){
            emiss_c_int[[i]][[q]][s,]	<- emiss_out$par
          }
        }
      }


      # Sample pouplaton values for gamma and conditional probabilities using Gibbs sampler -----------
      # gamma_mu0_int_n and gamma_mu_int_bar are matrices, with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
      gamma_mu0_int_n           <- solve(t(xx[[1]]) %*% xx[[1]] + gamma_K0_int)  %*% (t(xx[[1]]) %*% gamma_c_int[[i]] + gamma_K0_int %*% gamma_mu0_int[[i]])
      gamma_V_int_n             <- gamma_V0_int + t(gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_int_n) %*% (gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_int_n) + t(gamma_mu0_int_n - gamma_mu0_int[[i]]) %*% gamma_K0_int %*% (gamma_mu0_int_n - gamma_mu0_int[[i]])
      gamma_V_int[[i]]      <- solve(rwish(S = solve(gamma_V_int_n), v = gamma_nu_int + n_subj))
      gamma_mu_int_bar[[i]] <- gamma_mu0_int_n + solve(chol(t(xx[[1]]) %*% xx[[1]] + gamma_K0_int)) %*% matrix(rnorm((m - 1) * nx[1]), nrow = nx[1]) %*% t(solve(chol(solve(gamma_V_int[[i]]))))
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


      # Sample subject values for the intercepts of gamma using RW Metropolis sampler -----------
      for (s in 1:n_subj){
        gamma_int_candcov_comb 			<- chol2inv(chol(subj_data[[s]]$gamma_mhess[1:(m-1), 1:(m-1)] +
                                                     chol2inv(chol(gamma_V_int[[i]]))))
        gamma_int_RWout					    <- mnl_RW_once_int_fixed_bet(int1 = gamma_c_int[[i]][s,], Obs = c(trans[[s]][[i]], 1:m), n_cat = m,
                                                             xx_t =  matrix(c(rep(1,length(xx_t_state[[s]][[i]]) + 2), xx_t_state[[s]][[i]], rep(0,m)), ncol = 2),
                                                             bet = gamma_c_bet[[i]], mu_int_bar1 = c(t(gamma_mu_int_bar[[i]]) %*% xx[[1]][s,]),
                                                             V_int1 = gamma_V_int[[i]], scalar = gamma_scalar, candcov1 = gamma_int_candcov_comb)
        gamma_c_int[[i]][s,]		 <- gamma_int_RWout$draw_int
        gamma[[s]][i,] <- PD_subj[[s]][iter, c((sum(m * q_emiss) + 1 + (i - 1) * m):(sum(m * q_emiss) + (i - 1) * m + m))] <- gamma_int_RWout$prob
        gamma_int_naccept[s, i]			<- gamma_int_naccept[s, i] + gamma_int_RWout$accept
        gamma_int_subj[[s]][iter, (1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1))] <- gamma_c_int[[i]][s,]

        # Sample subject values for the conditional emission probabilities using RW Metropolis sampler -----------
        start <- c(0, q_emiss * m)
        for(q in 1:n_dep){
          emiss_candcov_comb		     <- chol2inv(chol(subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ] + chol2inv(chol(emiss_V_int[[i]][[q]]))))
          emiss_RWout				       <- mnl_RW_once(int1 = emiss_c_int[[i]][[q]][s,], Obs = c(cond_y[[s]][[i]][[q]], 1:q_emiss[q]), n_cat = q_emiss[q], mu_int_bar1 = c(t(emiss_mu_int_bar[[i]][[q]]) %*% xx[[1 + q]][s,]), V_int1 = emiss_V_int[[i]][[q]], scalar = emiss_scalar[[q]], candcov1 = emiss_candcov_comb)
          emiss[[s]][[q]][i,]		   <- PD_subj[[s]][iter, (sum(start[1:q]) + 1 + (i - 1) * q_emiss[q]):(sum(start[1:q]) + (i - 1) * q_emiss[q] + q_emiss[q])] <- emiss_RWout$prob
          emiss_naccept[[q]][s, i]	 <- emiss_naccept[[q]][s, i] + emiss_RWout$accept
          emiss_c_int[[i]][[q]][s,] <- emiss_RWout$draw_int
          emiss_int_subj[[s]][[q]][iter, (1 + (i - 1) * (q_emiss[q] - 1)) : ((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1))]	<- emiss_c_int[[i]][[q]][s,]
        }
      }
      # Sample fixed over subjects values for beta of gamma using RW Metropolis sampler -----------
      gamma_bet_candcov_comb 			<- chol2inv(chol(gamma_hess_mle_pooled[[i]][m:(2*(m-1)), m:(2*(m-1))] +
                                                   chol2inv(chol(gamma_V0_bet))))
      gamma_bet_RWout					    <- mnl_RW_once_fixed_bet_pooled(bet1 = gamma_c_bet[[i]], Obs_p = Obs_p, n_cat = m,
                                                           xx_t_p = xx_t_p, n_subj = n_subj, int_p = gamma_c_int[[i]],
                                                           mu_bet_bar1 = gamma_mu0_bet,
                                                           V_bet1 = gamma_V0_bet, scalar = gamma_scalar, candcov2 = gamma_bet_candcov_comb)
      gamma_c_bet[[i]]		<- gamma_bet_RWout$draw_int
      gamma_bet_naccept[i]			<- gamma_bet_naccept[i] + gamma_bet_RWout$accept

    }


    # End of 1 MCMC iteration, save output values --------
    gamma_int_bar[iter, ]				   	<- unlist(lapply(gamma_mu_int_bar, "[",1,))
    if(nx[1] > 1){
      gamma_cov_bar[iter, ]      	<- unlist(lapply(gamma_mu_int_bar, "[",-1,))
    }
    gamma_tv_cov[iter,]               <- unlist(gamma_c_bet)
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
    out <- list(input = list(m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,
                             burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                PD_subj = PD_subj, gamma_int_subj = gamma_int_subj, emiss_int_subj = emiss_int_subj,
                gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar, gamma_tv_cov = gamma_tv_cov,
                emiss_int_bar = emiss_int_bar, emiss_cov_bar = emiss_cov_bar, gamma_prob_bar = gamma_prob_bar,
                emiss_prob_bar = emiss_prob_bar, gamma_int_naccept = gamma_int_naccept, emiss_naccept = emiss_naccept,
                gamma_bet_naccept = gamma_bet_naccept, sample_path = sample_path)
  } else {
    out <- list(input = list(m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,
                             burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                PD_subj = PD_subj, gamma_int_subj = gamma_int_subj, emiss_int_subj = emiss_int_subj,
                gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar, gamma_tv_cov = gamma_tv_cov,
                emiss_int_bar = emiss_int_bar, emiss_cov_bar = emiss_cov_bar, gamma_prob_bar = gamma_prob_bar,
                emiss_prob_bar = emiss_prob_bar, gamma_int_naccept = gamma_int_naccept, emiss_naccept = emiss_naccept,
                gamma_bet_naccept = gamma_bet_naccept)
  }
  class(out) <- append(class(out), "mHMM_cat_tvc")
  return(out)
}
