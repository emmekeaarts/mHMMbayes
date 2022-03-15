#' Multilevel hidden  Markov model using Bayesian estimation for continuous
#' observations
#'
#' \code{mHMM_vary} fits a multilevel (also known as mixed or random effects)
#' hidden Markov model (HMM) to intense longitudinal data with categorical
#' and/or continuous (i.e., normally distributed) observations of multiple
#' subjects using Bayesian estimation, and creates an object of class mHMM_vary.
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
#' If covariates are specified and the user wants to set the values for the
#' parameters of the hyper-prior distributions manually, the specification of
#' the elements in the arguments of the hyper-prior parameter values of the
#' normal distribution on the means change as follows: the number of rows in the
#' matrices \code{gamma_mu0} and \code{emiss_mu0} are equal to 1 + the number of
#' covariates used to predict the transition probability matrix for
#' \code{gamma_mu0} and the emission distribution for \code{emiss_mu0} (i.e.,
#' the first row correspond to the hyper-prior mean values of the intercepts,
#' the subsequent rows correspond to the hyper-prior mean values of the
#' regression coefficients connected to each of the covariates), and
#' \code{gamma_K0} and \code{emiss_K0} are now a matrix with the number of
#' hypothetical prior subjects on which the vectors of means (i.e., the rows in
#' \code{gamma_mu0} or \code{emiss_mu0}) are based on the diagonal, and
#' off-diagonal elements equal to 0. Note that the hyper-prior parameter values
#' of the inverse Wishart distribution on the covariance matrix remains
#' unchanged, as the estimates of the regression coefficients for the covariates
#' are fixed over subjects.
#'
#'
#'
#'ADJUST: use of n_dep in the help file! change to n_cat and n_cont? introduce
#'
#' @param s_data A matrix containing the observations to be modelled, where the
#'   rows represent the observations over time. In \code{s_data}, the first
#'   column indicates subject id number. Hence, the id number is repeated over
#'   rows equal to the number of observations for that subject. The subsequent
#'   columns contain the dependent variable(s). Note that the dependent
#'   variables are assumed to be continuous (i.e., normally distributed
#'   depending on the hidden states). The total number of rows are equal to the
#'   sum over the number of observations of each subject, and the number of
#'   columns are equal to the number of dependent variables (\code{n_dep}) + 1.
#'   The number of observations can vary over subjects.
#' @param gen List containing the following elements denoting the general model
#'   properties: \itemize{
#'   \item{\code{m}: numeric vector with length 1 denoting the number of hidden
#'   states}
#'   \item{\code{n_dep}: numeric vector with length 1 denoting the number of
#'   dependent variables}
#'   \item{\code{q_emiss}: numeric vector with length \code{n_dep} denoting the
#'   number of observed categories for dependent variables with a categorical
#'   emission distribution. Input should equal 0 for dependent variables with a
#'   continous (i.e., normally distributed) emission distribution}}
#' @param data_distr String vector with either length 1 or length \code{n_dep}
#'   describing the observation type of the data to be simulated. Currently
#'   supported are \code{'categorical'} and \code{'continuous'}. When
#'   \code{'data_distr'} contains 1 value, the observation type is assumed equal
#'   over all dependent variables. The default equals to \code{data_distr =
#'   'categorical'}.
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
#'   values for gamma. The subsequent elements are matrices with \code{m} rows
#'   and 2 columns; the first column denoting the mean of state \emph{i} (row
#'   \emph{i}) and the second column denoting the variance of state \emph{i}
#'   (row \emph{i}) of the Normal distribution. Note that \code{start_val}
#'   should not contain nested lists (i.e., lists within lists).
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
#'  \itemize{\item{\code{gamma_mu0}: a list containing m matrices; one matrix
#'  for each row of the transition probability matrix gamma. Each matrix
#'  contains the hypothesized mean values of the intercepts. Hence, each matrix
#'  consists of one row (when not including covariates in the model) and
#'  \code{m} - 1 columns}
#'  \item{\code{gamma_K0}: numeric vector with length 1 denoting the number of
#'  hypothetical prior subjects on which the vector of means \code{gamma_mu0} is
#'  based}
#'  \item{\code{gamma_nu}: numeric vector with length 1 denoting the degrees of
#'  freedom of the Inverse Wishart distribution}
#'  \item{\code{gamma_V}: matrix of \code{m} - 1 by \code{m} - 1 containing the
#'  hypothesized variance-covariance matrix between the set of intercepts.}}
#'  Note that \code{gamma_K0}, \code{gamma_nu} and \code{gamma_V} are assumed
#'  equal over the states. The mean values of the intercepts (and regression
#'  coefficients of the covariates) denoted by \code{gamma_mu0} are allowed to
#'  vary over the states.
#'
#'  The default values for the hyper-prior on gamma are: all elements of the
#'  matrices contained in \code{gamma_mu0} set to 0, \code{gamma_K0} set to 1,
#'  \code{gamma_nu} set to 3 + m - 1, and the diagonal of \code{gamma_V} (i.e.,
#'  the variance) set to 3 + m - 1 and the off-diagonal elements (i.e., the
#'  covariance) set to 0.
#'
#'  See \emph{Details} below if covariates are used for changes in the settings
#'  of the arguments of \code{gamma_hyp_prior}.
#' @param emiss_cont_hyp_prior A list containing user specified parameters of
#'   the hyper-prior distribution on the continious dependent variables. That
#'   is, on the Normal (i.e., Gaussian) emission distributions
#'   (and regression coefficients given that covariates are used) for each of
#'   the states. The hyper-prior connected to the means of the Normal emission
#'   distribution(s) is a Normal-Inverse-Gamma distribution (i.e., assuming both
#'   unknown populaltion mean and variance between subject level means). The
#'   hyper-prior on each of fixed variances of the Normal emission distribuitons
#'   is an Inverse gamma distribution (i.e., assuming a known mean).
#'
#'  Hence, the list \code{emiss_hyp_prior} contains the following elements:
#'  \itemize{\item{\code{emiss_mu0}: a list containing \code{n_dep} matrices
#'  with one row (when not including covariates in the model) and \code{m}
#'  columns denoting the hypothesized mean values of the Normal emission
#'  distributions in each of the states for each dependent variable \code{k}}.
#'  \item{\code{emiss_K0}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables, where each element is an integer denoting
#'  the number of hypothetical prior subjects on which the vector of means
#'  \code{emiss_mu0} is based}.
#'  \item{\code{emiss_nu}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables, where each element is an integer
#'  denoting the degrees of freedom of the Inverse Gamma hyper-prior
#'  distribution connected to the emission distribution means (note: here, the
#'  Inverse Gamma hyper-prior distribution is parametrized as a scaled inverse
#'  chi-squared distribution).}
#'  \item{\code{emiss_V}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables \code{k}, where each element is a vector
#'  with lenght \code{m} containing the hypothesized variances between the
#'  between subject level means of the Inverse Gamma hyper-prior distribution
#'  connected to the emission distribution means (note: here, the Inverse Gamma
#'  hyper-prior distribution is parametrized as a scaled inverse chi-squared
#'  distribution).}
#'  \item{\code{emiss_a0}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables \code{k}, where each element is a vector
#'  with lenght \code{m} containing the shape values of the Inverse Gamma
#'  hyper-prior on each of fixed variances of the Normal emission distribuitons
#'  (note: here the standard Inverse Gamma parametrization is used).}
#'  \item{\code{emiss_b0}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables \code{k}, where each element is a vector
#'  with lenght \code{m} containing the scale values of the Inverse Gamma
#'  hyper-prior on each of fixed variances of the Normal emission distribuitons
#'  (note: here the standard Inverse Gamma parametrization is used).}}
#'  Note that \code{emiss_K0} and \code{emiss_nu} are assumed
#'  equal over the states.
#'
#'
#'  See \emph{Details} below if covariates are used for changes in the settings
#'  of the arguments of \code{emiss_hyp_prior}.
#'
#'  @param emiss_cat_hyp_prior An optional list containing user specified
#'    parameters of the hyper prior distribution on the categorical dependent
#'    variables. That is, on the multivariate normal distribution of the
#'    intercepts (and regression coefficients given that covariates are used) of
#'    the multinomial regression model of the emission distribution. The
#'    hyper-prior for the mean intercepts is a multivariate Normal distribution,
#'    the hyper-prior for the covariance matrix between the set of (state
#'    specific) intercepts is an Inverse Wishart distribution.
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
#'  and the diagonal of \code{gamma_V} (i.e., the variance) set to 3 +
#'  \code{q_emiss[k]} - 1 and the off-diagonal elements (i.e., the covariance)
#'  set to 0.
#'
#'  See \emph{Details} below if covariates are used for changes in the settings
#'  of the arguments of \code{emiss_hyp_prior}.
#' @param emiss_sampler An optional list containing user specified settings for
#'   the proposal distribution of the random walk (RW) Metropolis sampler for
#'   the subject level parameter estimates of the intercepts modeling the
#'   emission distributions of the categorical dependent variables \code{k}. The list
#'   \code{emiss_sampler} contains the following elements:
#'  \itemize{\item{\code{emiss_int_mle0}: a list containing \code{n_dep} elements
#'  corresponding to each of the categorical dependent variables \code{k}, where each
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
#'  @param gamma_sampler An optional list containing user specified settings for
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
#'
#' @return \code{mHMM_cont} returns an object of class \code{mHMM_cont}, which has
#'   \code{print} and \code{summary} methods to see the results.
#'   The object contains the following components:
#'   \describe{
#'   \item{\code{PD_subj}}{A list containing one list per subject with the
#'   elements \code{trans_prob}, \code{cat_emiss}, \code{cont_emiss}, and
#'   \code{log_likl} providing the subject parameter esimates over the
#'   iterations of the MCMC sampler of the transition probabilities gamma, the
#'   categorical emission distribution (emission probabilities), and the
#'   continious emission distribuituion (subsequently the the emission means and
#'   the (fixed over subjects) emission variances), respectiviely, and the log
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
#'   \item{\code{gamma_naccept}}{A matrix containing the number of accepted
#'   draws at the subject level RW Metropolis step for each set of parameters of
#'   the transition probabilities. The subjects are contained in the rows, and
#'   the columns contain the sets of parameters.}
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
#'   Normal emision distributions. over the iterations of the Gibbs sampler. The
#'   iterations of the sampler are contained in the rows of the matrix, and the
#'   columns contain the group level variance in the mean.}
#'   \item{\code{emiss_var_bar}}{A list containing one matrix per dependent
#'   variable, denoting the (fixed over subjects) variance of the Normal emision
#'   distributions over the iterations of the Gibbs sampler. The iterations of
#'   the sampler are contained in the rows of the matrix, and the columns
#'   contain the group level emission variances.}
#'   \item{\code{emiss_cont_cov_bar}}{A list containing one matrix per
#'   continuous dependent variable, denoting the group level regression
#'   coefficients predicting the emission means within each of the continuous
#'   dependent variables over the iterations of the Gibbs sampler. The
#'   iterations of the sampler are contained in the rows  of the matrix, and the
#'   columns contain the group level regression coefficients.}
#'   \item{\code{emiss_cat_cov_bar}}{A list containing one matrix per
#'   categorical dependent variable, denoting the group level regression
#'   coefficients of the multinomial logistic regression predicting the emission
#'   probabilities within each of the categorical dependent variables over the
#'   iterations of the Gibbs sampler. The iterations of the sampler are
#'   contained in the rows  of the matrix, and the columns contain the group
#'   level regression coefficients.}
#'   \item{\code{label_switch}}{A matrix of \code{m} rows and \code{n_dep}
#'   columns containing the percentage of times the group mean of the emission
#'   distriubion of state \code{i} was sampled to be a smaller value compared to
#'   the group mean of of the emission distriubion of state \code{i-1}. If the
#'   state dependent means of the emission distributions were given in a ranked
#'   order (low to high) to both the start values and hyper-priors, a high
#'   percentage in \code{label_switch} indicates that label switching possibly
#'   poses a problem in the analysis, and further diagnostics (e.g.,
#'   traceplots and posterior distributions) should be inspected.}
#'   \item{\code{input}}{Overview of used input specifications: the number of
#'   states \code{m}, the number of used dependent variables \code{n_dep}, the
#'   number of iterations \code{J} and the specified burn in period
#'   \code{burn_in} of the hybrid Metropolis within Gibbs sampler, the number of
#'   subjects \code{n_subj}, the observation length for each subject
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
#' ###### Example on simulated data
#' # simulating data with both categorical and continous dependent variables
#'
#' number of observations per subject, subjects, states and dependent variables
#' n_t     <- 100
#' n       <- 10
#' m       <- 3
#' n_dep   <- 3
#'
#' # distribution of each of the dependent variables
#' data_distr = c("categorical", "continuous", "continuous")
#'
#' # number of categories for the categorical depenent variable
#' q_emiss <- c(4, 0, 0)
#'
#' # transition probabilities between the states
#' gamma <- matrix(c(0.6, 0.3, 0.1,
#'                   0.2, 0.7, 0.1,
#'                   0.1, 0.4, 0.5), nrow = m, byrow = TRUE)
#'
#' # emission distributions
#'                     # observation probabilities of each of the 4 categories of the first (categorical) dependent variable over the 3 states
#' emiss_distr <- list(matrix(c(0.80, 0.13, 0.05, 0.02,
#'                     0.08, 0.85, 0.05, 0.02,
#'                     0.05, 0.20, 0.40, 0.35), nrow = m, byrow = TRUE),
#'                     # mean and variance of second (continuous) dependend in each of the 3 states
#'                     matrix(c(6.0, 0.5,
#'                              5.0, 0.5,
#'                              3.0, 0.5), nrow = m, byrow = TRUE),
#'                     # mean and variance of third (continous) in each of the 3 states
#'                     matrix(c(4.0, 1.0,
#'                            5.0, 1.0,
#'                            6.5, 1.0), nrow = m, byrow = TRUE))
#'
#' # amount of variance in the transition probabilities between subjects
#' var_gamma <- 0.01
#'
#' # amount of variance in state dependent distributions of the three dependent variables between subjects
#' var_emiss <- c(0.01, 0.10, 0.10)
#'
#' # simulate the data
#' data_vary <- sim_mHMM(n_t = n_t, n = n, data_distr = data_distr, m = m , n_dep = n_dep, q_emiss = q_emiss,
#'                   gamma = gamma, emiss_distr = emiss_distr, var_gamma = var_gamma, var_emiss = var_emiss)
#'
#' # Specify hyper-prior for the continuous emission distribution
#' hyp_pr = emiss_cont_hyp_prior = list(
#'            emiss_mu0 = list(matrix(c(6,5,3), nrow = 1),
#'                             matrix(c(4,5,6), nrow = 1)),
#'            emiss_K0 = list(1,1),
#'            emiss_nu = list(1,1),
#'            emiss_V = list(rep(2, 3), rep(2, 3)),
#'            emiss_a0 = list(rep(1, 3), rep(1, 3)),
#'            emiss_b0 = list(rep(1, 3), rep(1, 3))
#'            )
#'
#' # Run the model on the simulated data:
#' out_3st_vary_dep <- mHMM_vary(s_data = data_vary$obs,
#'                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                     data_distr = data_distr,
#'                     start_val = c(list(gamma), emiss_distr),
#'                     emiss_cont_hyp_prior = hyp_pr,
#'                     mcmc = list(J = 11, burn_in = 5))
#'
#'
#' @export
#'
#'

mHMM_vary <- function(s_data, gen, data_distr, xx = NULL, start_val, emiss_cont_hyp_prior, mcmc, return_path = FALSE, print_iter, show_progress = TRUE,
                      emiss_cat_hyp_prior = NULL, emiss_sampler = NULL, gamma_hyp_prior = NULL, gamma_sampler = NULL){

  if(!missing(print_iter)){
    warning("The argument print_iter is depricated; please use show_progress instead to show the progress of the algorithm.")
  }
  # Initialize data -----------------------------------
  # dependent variable(s), sample size, dimensions gamma and conditional distribuiton
  n_dep			 <- gen$n_dep
  n_cat      <- sum(data_distr %in% 'categorical')
  n_cont     <- sum(data_distr %in% 'continuous')
  which_cat  <- which(data_distr %in% 'categorical')
  which_cont <- which(data_distr %in% 'continuous')
  # include check that lenght data_distr is either one, or equal to n_dep
  if(n_dep > 1 & length(data_distr) == 1){
    data_distr <- rep(data_distr, n_dep)
  }
  dep_labels <- colnames(s_data[,2:(n_dep+1)])
  id         <- unique(s_data[,1])
  n_subj     <- length(id)
  subj_data  <- rep(list(NULL), n_subj)
  if(sum(sapply(s_data, is.factor)) > 0 ){
    stop("Your data contains factorial variables, which cannot be used as input in the function mHMM. All variables have to be numerical.")
  }
  for(s in 1:n_subj){
    subj_data[[s]]$y <- as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep)
  }
  ypooled    <- n <- NULL
  n_vary     <- numeric(n_subj)
  m          <- gen$m
  if(n_cat > 0){
    q_emiss 		 <- gen$q_emiss
    emiss_int_mle <- rep(list(NULL), n_dep)
    emiss_mhess   <- rep(list(NULL), n_dep)
    for(q in which_cat){
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
  } else {
    for(s in 1:n_subj){
      ypooled   <- rbind(ypooled, subj_data[[s]]$y)
      n         <- dim(subj_data[[s]]$y)[1]
      n_vary[s] <- n
      subj_data[[s]]	<- c(subj_data[[s]], n = n, list(gamma_converge = numeric(m), gamma_int_mle = matrix(, m, (m - 1)),
                                                      gamma_mhess = matrix(, (m - 1) * m, (m - 1))))
    }
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
    gamma_scalar    <- 2.93 / sqrt(m - 1)
    gamma_w         <- .1
  } else {
    gamma_int_mle0  <- gamma_sampler$gamma_int_mle0
    gamma_scalar    <- gamma_sampler$gamma_scalar
    gamma_w         <- gamma_sampler$gamma_w
  }



  # Initialize Gamma hyper prior
  if(is.null(gamma_hyp_prior)){
    gamma_mu0	  <- rep(list(matrix(0,nrow = nx[1], ncol = m - 1)), m)
    gamma_K0			<- diag(1, nx[1])
    gamma_nu			<- 3 + m - 1
    gamma_V			  <- gamma_nu * diag(m - 1)
  } else {
    ###### BUILD in a warning / check if gamma_mu0 is a matrix when given, with  nrows equal to the number of covariates
    gamma_mu0			<- gamma_hyp_prior$gamma_mu0
    gamma_K0			<- diag(gamma_hyp_prior$gamma_K0, nx[1])
    gamma_nu			<- gamma_hyp_prior$gamma_nu
    gamma_V			  <- gamma_hyp_prior$gamma_V
  }


  # Initialize continous emiss hyper prior
  if(n_cont > 0){
    if(missing(emiss_cont_hyp_prior)){
      stop("The hyper-prior values for the Normal emission distribution(s) denoted by emiss_cont_hyp_prior needs to be specified")
    }

    # emiss_mu0: a list containing n_cont matrices with in the first row the hypothesized mean values of the Normal emission
    # distributions in each of the states over the m coloumns. Subsequent rows contain the hypothesised regression
    # coefficients for covariates influencing the state dependent mean value of the normal distribution
    emiss_cont_mu0	  <- rep(list(NULL), n_cont)
    emiss_cont_a0	  <- rep(list(NULL), n_cont)
    emiss_cont_b0	  <- rep(list(NULL), n_cont)
    emiss_cont_V	  <- rep(list(NULL), n_cont)
    emiss_cont_nu	    <- rep(list(NULL), n_cont)
    emiss_cont_K0     <- rep(list(NULL), n_cont)
    for(q in 1:n_cont){
      # emiss_hyp_prior[[q]]$emiss_mu0 has to contain a list with lenght equal to m, and each list contains matrix with number of rows equal to number of covariates for that dep. var.
      # stil build in a CHECK, with warning / stop / switch to default prior
      emiss_cont_mu0[[q]]	 <- emiss_cont_hyp_prior$emiss_mu0[[q]]
      emiss_cont_nu[[q]]	 <- emiss_cont_hyp_prior$emiss_nu[[q]]
      emiss_cont_V[[q]]		 <- emiss_cont_hyp_prior$emiss_V[[q]]
      emiss_cont_K0[[q]]	 <- diag(emiss_cont_hyp_prior$emiss_K0, nx[1 + q])
      emiss_cont_a0[[q]] <- emiss_cont_hyp_prior$emiss_a0[[q]]
      emiss_cont_b0[[q]] <- emiss_cont_hyp_prior$emiss_b0[[q]]
    }
  }

  # Initialize categorical emiss sampler
  if(n_cat > 0){
    if(is.null(emiss_sampler)){
      emiss_int_mle0 <- rep(list(NULL), n_cat)
      emiss_scalar	<- rep(list(NULL), n_cat)
      for(q in 1:n_cat){
        emiss_int_mle0[[q]] <- rep(0, q_emiss[which_cat][q] - 1)
        emiss_scalar[[q]] 	<- 2.93 / sqrt(q_emiss[which_cat][q] - 1)
      }
      emiss_w		<- .1
    } else {
      emiss_int_mle0	<- emiss_sampler$emiss_int_mle0
      emiss_scalar 	<- emiss_sampler$emiss_scalar
      emiss_w    		<- emiss_sampler$emiss_w
    }

    # Initialize Pr hyper prior
    # emiss_mu0: for each dependent variable, emiss_mu0 is a list, with one element for each state.
    # Each element is a matrix, with number of rows equal to the number of covariates (with the intercept being one cov),
    # and the number of columns equal to q_emiss[q] - 1.
    emiss_cat_mu0	  <- rep(list(vector("list", m)), n_cat)
    emiss_cat_nu	    <- rep(list(NULL), n_cat)
    emiss_cat_V	    <- rep(list(NULL), n_cat)
    emiss_cat_K0     <- rep(list(NULL), n_cat)
    if(is.null(emiss_cat_hyp_prior)){
      for(q in 1:n_cat){
        ind <- which_cat[q]
        for(i in 1:m){
          emiss_cat_mu0[[q]][[i]]		<- matrix(0, ncol = q_emiss[ind] - 1, nrow = nx[1 + ind])
        }
        emiss_cat_nu[[q]]		<- 3 + q_emiss[ind] - 1
        emiss_cat_V[[q]]		  <- emiss_cat_nu[[q]] * diag(q_emiss[ind] - 1)
        emiss_cat_K0[[q]]		<- diag(1, nx[1 + ind])
      }
    } else {
      for(q in 1:n_cat){
        # emiss_hyp_prior[[q]]$emiss_mu0 has to contain a list with lenght equal to m, and each list contains matrix with number of rows equal to number of covariates for that dep. var.
        # stil build in a CHECK, with warning / stop / switch to default prior
        emiss_cat_mu0[[q]]	 <- emiss_cat_hyp_prior$emiss_mu0[[q]]
        emiss_cat_nu[[q]]	 <- emiss_cat_hyp_prior$emiss_nu[[q]]
        emiss_cat_V[[q]]		 <- emiss_cat_hyp_prior$emiss_V[[q]]
        emiss_cat_K0[[q]]	 <- diag(emiss_cat_hyp_prior$emiss_K0, nx[1 + which_cat[q]])
      }
    }
  }


#####################

  # Define objects used to store data in mcmc algorithm, not returned ----------------------------
  # overall
  c <- llk <- numeric(1)
  sample_path <- lapply(n_vary, dif_matrix, cols = J)
  trans <- rep(list(vector("list", m)), n_subj)

  # gamma
  gamma_int_mle_pooled <- gamma_pooled_ll <- vector("list", m)
  gamma_c_int <- rep(list(matrix(, n_subj, (m-1))), m)
  gamma_mu_int_bar <- gamma_V_int <- vector("list", m)
  gamma_mu_prob_bar <- rep(list(numeric(m)), m)
  gamma_naccept <- matrix(0, n_subj, m)

  # emiss
  cond_y <- lapply(rep(n_dep, n_subj), nested_list, m = m)
  cond_y_pooled <- rep(list(rep(list(NULL),n_dep)), m)

    # categorical
    emiss_int_mle_pooled <- emiss_pooled_ll <- rep(list(vector("list", n_cat)), m)
    emiss_c_int <- rep(list(lapply(q_emiss[which_cat] - 1, dif_matrix, rows = n_subj)), m)
    emiss_mu_int_bar <- emiss_V_int <- rep(list(vector("list", n_cat)), m)
    emiss_mu_prob_bar <- rep(list(lapply(q_emiss[which_cat], dif_vector)), m)
    emiss_naccept <- rep(list(matrix(0, n_subj, m)), n_cat)

    # continuous
    emiss_c_mu <- rep(list(rep(list(matrix(,ncol = 1, nrow = n_subj)),n_cont)), m)
    for(i in 1:m){
     for(q in 1:n_cont){
       ind <- which_cont[q]
       emiss_c_mu[[i]][[q]][,1] <- start_val[[1 + ind]][i,1]
     }
    }
    emiss_V_mu <- emiss_c_mu_bar <- emiss_c_V <- rep(list(rep(list(NULL),n_cont)), m)
    ss_subj <- n_cond_y <- numeric(n_subj)
    label_switch <- matrix(0, ncol = n_cont, nrow = m, dimnames = list(c(paste("mu_S", 1:m, sep = "")), dep_labels[which_cont]))


  # Define objects that are returned from mcmc algorithm ----------------------------
  # Define object for subject specific posterior density, put start values on first rows
  if(length(start_val) != n_dep + 1){
    stop("The number of elements in the list start_val should be equal to 1 + the number of dependent variables,
         and should not contain nested lists (i.e., lists within lists)")
  }

  PD              <- list(trans_prob = matrix(, nrow = J, ncol = m * m),
                          cat_emiss = if(n_cat > 0){
                            matrix(, nrow = J, ncol = sum(m * q_emiss))} else {
                              NULL
                            },
                          cont_emiss = if(n_cont) {
                            matrix(, nrow = J, ncol = (n_cont * 2 * m))} else {
                              NULL
                            },
                          log_likl = matrix(, nrow = J, ncol = 1))
  colnames(PD$trans_prob) 	<- paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = "")
  PD$trans_prob[1, ] <- unlist(sapply(start_val, t))[1:(m*m)]
  if (n_cat > 0){
    PD_cat_emiss_names   <- character()
    for(q in which_cat){
      PD_cat_emiss_names <- c(PD_cat_emiss_names, paste("dep", q, "_S", rep(1:m, each = q_emiss[q]), "_emiss", rep(1:q_emiss[q], m), sep = ""))
    }
    colnames(PD$cat_emiss)    <- PD_cat_emiss_names
    for(q in 1:n_cat){
      ind <- which_cat[q]
      PD$cat_emiss[1, (sum(c(0,q_emiss)[1 : ind] * m) + 1):sum(q_emiss[1 : ind] * m)] <- as.vector(t(start_val[[ind + 1]]))
    }
  }
  if(n_cont > 0){
    colnames(PD$cont_emiss) <- c(paste("dep", rep(which_cont, each = m), "_S", rep(1:m), "_mu", sep = ""),
                                 paste("dep", rep(which_cont, each = m), "_S", rep(1:m), "_fixvar", sep = ""))
    for(q in 1:n_cont){
      ind <- which_cont[q]
      PD$cont_emiss[1, ((q-1) * m + 1):(q * m)] <- start_val[[ind + 1]][,1]
      PD$cont_emiss[1,  (n_cont * m + (q-1) * m + 1):(n_cont * m + q * m)] <- start_val[[ind + 1]][,2]
    }
  }
  colnames(PD$log_likl) <-  "LL"
  PD_subj				<- rep(list(PD), n_subj)

  # Define object for population posterior density (probabilities and regression coefficients parameterization )
  gamma_prob_bar		<- matrix(, nrow = J, ncol = (m * m))
  colnames(gamma_prob_bar) <- paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = "")
  gamma_prob_bar[1,] <- PD$trans_prob[1,]
  gamma_int_bar				<- matrix(, nrow = J, ncol = ((m-1) * m))
  colnames(gamma_int_bar) <- paste("int_S", rep(1:m, each = m-1), "toS", rep(2:m, m), sep = "")
  gamma_int_bar[1,] <- as.vector(prob_to_int(matrix(gamma_prob_bar[1,], byrow = TRUE, ncol = m, nrow = m)))
  if(nx[1] > 1){
    gamma_cov_bar				<- matrix(, nrow = J, ncol = ((m-1) * m) * (nx[1] - 1))
    # colnames(gamma_cov_bar) <- paste( paste("cov", rep(1 : (nx[1] - 1),each = nx[1]-1), "_", sep = ""), "S", rep(1:m, each = (m-1) * (nx[1] - 1)), "toS", rep(2:m, m * (nx[1] - 1)), sep = "")
    colnames(gamma_cov_bar) <- paste( paste("cov", rep(1 : (nx[1] - 1),each = (m-1)), "_", sep = ""), "S", rep(1:m, each = (m-1) * (nx[1] - 1)), "toS", rep(2:m, m * (nx[1] - 1)), sep = "")
    gamma_cov_bar[1,] <- 0
  } else{
    gamma_cov_bar <- "No covariates where used to predict the transition probability matrix"
  }
  # Define object for subject specific posterior density (regression coefficients parameterization )
  gamma_int_subj			<- rep(list(gamma_int_bar), n_subj)

# continuous
  if(n_cont > 0){
    emiss_mu_bar			<- rep(list(matrix(, ncol = m, nrow = J, dimnames = list(NULL, c(paste("mu_", 1:m, sep = ""))))), n_cont)
    names(emiss_mu_bar) <- dep_labels[which_cont]
    for(q in 1:n_cont){
      emiss_mu_bar[[q]][1,] <- PD$cont_emiss[1, ((q-1) * m + 1):(q * m)]
    }
    emiss_varmu_bar			<- rep(list(matrix(, ncol = m, nrow = J, dimnames = list(NULL, c(paste("varmu_", 1:m, sep = ""))))), n_cont)
    names(emiss_varmu_bar) <- dep_labels[which_cont]
    emiss_var_bar			<- rep(list(matrix(, ncol = m, nrow = J, dimnames = list(NULL, c(paste("var_", 1:m, sep = ""))))), n_cont)
    names(emiss_var_bar) <- dep_labels[which_cont]
    for(q in 1:n_cont){
      emiss_var_bar[[q]][1,] <- PD$cont_emiss[1, (n_cont * m + (q-1) * m + 1):(n_cont * m + q * m)]
    }
    if(sum(c(nx[-1] - 1) * c(data_distr %in% 'continuous')) > 0){
      emiss_cont_cov_bar			<- lapply(m * (nx[-1][which_cont] - 1 ), dif_matrix, rows = J)
      names(emiss_cont_cov_bar) <- dep_labels
      for(q in 1:n_cont){
        ind <- which_cont[q]
        if(nx[1 + ind] > 1){
          # colnames(emiss_cont_cov_bar[[q]]) <-  paste( paste("cov", rep(1 : (nx[1+q] - 1),each = nx[1+q]-1), "_", sep = ""), "mu_S", rep(1:m, each = (nx[1 + q] - 1)), sep = "")
          colnames(emiss_cont_cov_bar[[q]]) <-  paste( paste("cov", 1 : (nx[1+ind] - 1), "_", sep = ""), "mu_S", rep(1:m, each = (nx[1 + ind] - 1)), sep = "")
          emiss_cont_cov_bar[[q]][1,] <- 0
        } else {
          emiss_cont_cov_bar[[q]] <- "No covariates where used to predict the emission probabilities for this outcome"
        }
      }
    } else{
      emiss_cont_cov_bar <- "No covariates where used to predict the continious emission probabilities"
    }
  }

  # categorical
  if(n_cat > 0){
    emiss_prob_bar			<- lapply(q_emiss[which_cat] * m, dif_matrix, rows = J)
    names(emiss_prob_bar) <- dep_labels[which_cat]
    for(q in 1:n_cat){
      ind <- which_cat[q]
      colnames(emiss_prob_bar[[q]]) <- paste("S", rep(1:m, each = q_emiss[ind]), "_emiss", rep(1:q_emiss[ind], m), sep = "")
      start <- c(0, q_emiss[which_cat] * m)
      emiss_prob_bar[[q]][1,] <- PD$cat_emiss[1,(sum(start[1:q]) + 1):(sum(start[1:q]) + (m * q_emiss[q]))]
    }
    emiss_int_bar			<- lapply((q_emiss[which_cat]-1) * m, dif_matrix, rows = J)
    names(emiss_int_bar) <- dep_labels[which_cat]
    for(q in 1:n_cat){
      ind <- which_cat[q]
      colnames(emiss_int_bar[[q]]) <-  paste("S", rep(1:m, each = q_emiss[ind] - 1), "_int_emiss", rep(2:q_emiss[ind], m), sep = "")
      emiss_int_bar[[q]][1,] <- as.vector(prob_to_int(matrix(emiss_prob_bar[[q]][1,], byrow = TRUE, ncol = q_emiss[ind], nrow = m)))
    }
    if(sum(c(nx[-1] - 1) * c(data_distr %in% 'categorical')) > 0){
      emiss_cat_cov_bar			<- lapply((q_emiss[which_cat]-1) * m * (nx[-1][which_cat] - 1 ), dif_matrix, rows = J)
      names(emiss_cat_cov_bar) <- dep_labels[which_cat]
      for(q in 1:n_cat){
        ind <- which_cat[q]
        if(nx[1 + ind] > 1){
          colnames(emiss_cat_cov_bar[[q]]) <-  paste( paste("cov", 1 : (nx[1 + ind] - 1), "_", sep = ""), "_S", rep(1:m, each = (q_emiss[ind] - 1) * (nx[1 + ind] - 1)), "_emiss", rep(2:q_emiss[ind], m * (nx[1 + ind] - 1)), sep = "")
          emiss_cat_cov_bar[[q]][1,] <- 0
        } else {
          emiss_cat_cov_bar[[q]] <- "No covariates where used to predict the emission probabilities for this outcome"
        }
      }
    } else{
      emiss_cat_cov_bar <- "No covariates where used to predict the categorical emission probabilities"
    }
    # Define object for subject specific posterior density (regression coefficients parameterization )
    emiss_int_subj			<- rep(list(emiss_int_bar), n_subj)
  }

  # Put starting values in place for fist run forward algorithm

  emiss				<- rep(list(start_val[2:(n_dep + 1)]), n_subj)
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
      # Run forward algorithm, obtain subject specific forward proababilities and log likelihood
      forward				<- vary_mult_fw_r_to_cpp(x = subj_data[[s]]$y, m = m, emiss = emiss[[s]],
                                          gamma = gamma[[s]], n_dep = n_dep, delta=NULL, data_distr = data_distr)
      alpha         <- forward[[1]]
      c             <- max(forward[[2]][, subj_data[[s]]$n])
      llk           <- c + log(sum(exp(forward[[2]][, subj_data[[s]]$n] - c)))
      PD_subj[[s]]$log_likl[iter, 1] <- llk

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
          cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q])
        }
      }
    }

    # The remainder of the mcmc algorithm is state specific
    for(i in 1:m){

      # Obtain MLE of the covariance matrices and log likelihood of gamma and emiss at subject and population level -----------------
      # used to scale the propasal distribution of the RW Metropolis sampler

      # population level, transition matrix
      trans_pooled			  <- factor(c(unlist(sapply(trans, "[[", i)), c(1:m)))
      gamma_mle_pooled		<- optim(gamma_int_mle0, llmnl_int, Obs = trans_pooled,
                                 n_cat = m, method = "BFGS", hessian = TRUE,
                                 control = list(fnscale = -1))
      gamma_int_mle_pooled[[i]]  <- gamma_mle_pooled$par
      gamma_pooled_ll[[i]]			<- gamma_mle_pooled$value

      # subject level, transition matrix
      for (s in 1:n_subj){
        wgt 				<- subj_data[[s]]$n / n_total

        # subject level, transition matrix
        gamma_out					<- optim(gamma_int_mle_pooled[[i]], llmnl_int_frac, Obs = c(trans[[s]][[i]], c(1:m)),
                               n_cat = m, pooled_likel = gamma_pooled_ll[[i]], w = gamma_w, wgt = wgt,
                               method="BFGS", hessian = TRUE, control = list(fnscale = -1))
        if(gamma_out$convergence == 0){
          subj_data[[s]]$gamma_converge[i] <- 1
          subj_data[[s]]$gamma_int_mle[i,] <- gamma_out$par
          subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ]	<-
            mnlHess_int(int = gamma_out$par, Obs = c(trans[[s]][[i]], c(1:m)), n_cat =  m)
        } else {
          subj_data[[s]]$gamma_converge[i] <- 0
          subj_data[[s]]$gamma_int_mle[i,] <- rep(0, m - 1)
          subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ]	<- diag(m-1)
        }

        # if this is first iteration, use MLE for current values RW metropolis sampler
        if (iter == 2){
          gamma_c_int[[i]][s,]		<- gamma_out$par
        }
      }

      # emission probabilities, seperate for each dependent variable, only for CATEGORICAL dependent variables
      if(n_cat > 0){
        for(q in 1:n_cat){
          ind <- which_cat[q]

          # population level
          cond_y_pooled					      <- numeric()
          ### MOET OOK ECHT BETER KUNNEN, eerst # cond_y_pooled				<- unlist(sapply(cond_y, "[[", m))
          for(s in 1:n_subj){
            cond_y_pooled             <- c(cond_y_pooled, cond_y[[s]][[i]][[ind]])
          }
          emiss_mle_pooled		<- optim(emiss_int_mle0[[q]], llmnl_int, Obs = c(cond_y_pooled, c(1:q_emiss[ind])),
                                     n_cat = q_emiss[ind], method = "BFGS", hessian = TRUE,
                                     control = list(fnscale = -1))
          emiss_int_mle_pooled[[i]][[q]]  <- emiss_mle_pooled$par
          emiss_pooled_ll[[i]][[q]]				<- emiss_mle_pooled$value

          # subject level
          for (s in 1:n_subj){
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
      }

      # Sample pouplaton values for gamma and conditional probabilities using Gibbs sampler -----------
      # gamma_mu0_n and gamma_mu_int_bar are matrices, with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
      gamma_mu0_n           <- solve(t(xx[[1]]) %*% xx[[1]] + gamma_K0)  %*% (t(xx[[1]]) %*% gamma_c_int[[i]] + gamma_K0 %*% gamma_mu0[[i]])
      gamma_V_n             <- gamma_V + t(gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) %*% (gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) + t(gamma_mu0_n - gamma_mu0[[i]]) %*% gamma_K0 %*% (gamma_mu0_n - gamma_mu0[[i]])
      gamma_V_int[[i]]      <- solve(rwish(S = solve(gamma_V_n), v = gamma_nu + n_subj))
      gamma_mu_int_bar[[i]] <- gamma_mu0_n + solve(chol(t(xx[[1]]) %*% xx[[1]] + gamma_K0)) %*% matrix(rnorm((m - 1) * nx[1]), nrow = nx[1]) %*% t(solve(chol(solve(gamma_V_int[[i]]))))
      gamma_exp_int				  <- matrix(exp(c(0, gamma_mu_int_bar[[i]][1,] )), nrow  = 1)
      gamma_mu_prob_bar[[i]] 	<- gamma_exp_int / as.vector(gamma_exp_int %*% c(rep(1,(m))))

      # Sample subject values for transition matrix -----------
      for (s in 1:n_subj){
        # Sample subject values for gamma using RW Metropolis sampler   ---------
        gamma_candcov_comb 			<- chol2inv(chol(subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ] + chol2inv(chol(gamma_V_int[[i]]))))
        gamma_RWout					    <- mnl_RW_once(int1 = gamma_c_int[[i]][s,], Obs = trans[[s]][[i]], n_cat = m, mu_int_bar1 = c(t(gamma_mu_int_bar[[i]]) %*% xx[[1]][s,]), V_int1 = gamma_V_int[[i]], scalar = gamma_scalar, candcov1 = gamma_candcov_comb)
        gamma[[s]][i,]  	      <- PD_subj[[s]]$trans_prob[iter, ((i-1) * m + 1) : ((i-1) * m + m)] <- gamma_RWout$prob
        gamma_naccept[s, i]			<- gamma_naccept[s, i] + gamma_RWout$accept
        gamma_c_int[[i]][s,]		<- gamma_RWout$draw_int
        gamma_int_subj[[s]][iter, (1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1))] <- gamma_c_int[[i]][s,]

        if(i == m){
          delta[[s]] 		<- solve(t(diag(m) - gamma[[s]] + 1), rep(1, m))
        }
      }

      # sample parameters for the emission distributions
      # note: the posterior is thus one of a Bayesian linear regression because of the optional regression parameters
      if(n_cat > 0){
        # categorical dependent variables:
        start <- c(0, q_emiss[which_cat] * m)

        for(q in 1:n_cat){
          ind <- which_cat[q]

          # population parameters
          emiss_mu0_n                 <- solve(t(xx[[1 + ind]]) %*% xx[[1 + ind]] + emiss_cat_K0[[q]]) %*% (t(xx[[1 + ind]]) %*% emiss_c_int[[i]][[q]] + emiss_cat_K0[[q]] %*% emiss_cat_mu0[[q]][[i]])
          emiss_V_n                   <- emiss_cat_V[[q]] + t(emiss_c_int[[i]][[q]] - xx[[1 + ind]] %*% emiss_mu0_n) %*% (emiss_c_int[[i]][[q]] - xx[[1 + ind]] %*% emiss_mu0_n) + t(emiss_mu0_n - emiss_cat_mu0[[q]][[i]]) %*% emiss_cat_K0[[q]] %*% (emiss_mu0_n - emiss_cat_mu0[[q]][[i]])
          emiss_V_int[[i]][[q]]       <- solve(rwish(S = solve(emiss_V_n), v = emiss_cat_nu[[q]] + n_subj))
          emiss_mu_int_bar[[i]][[q]]	 <- emiss_mu0_n + solve(chol(t(xx[[1 + ind]]) %*% xx[[1 + ind]] + emiss_cat_K0[[q]])) %*% matrix(rnorm((q_emiss[ind] - 1) * nx[1 + ind]), nrow = nx[1 + ind]) %*% t(solve(chol(solve(emiss_V_int[[i]][[q]]))))
          emiss_exp_int				       <- matrix(exp(c(0, emiss_mu_int_bar[[i]][[q]][1, ])), nrow  = 1)
          emiss_mu_prob_bar[[i]][[q]] <- emiss_exp_int / as.vector(emiss_exp_int %*% c(rep(1, (q_emiss[ind]))))

          # subject level parameters sampled using RW Metropolis sampler
          for (s in 1:n_subj){
            emiss_candcov_comb		     <- chol2inv(chol(subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[ind] - 1)):((q_emiss[ind] - 1) + (i - 1) * (q_emiss[ind] - 1)), ] + chol2inv(chol(emiss_V_int[[i]][[q]]))))
            emiss_RWout				       <- mnl_RW_once(int1 = emiss_c_int[[i]][[q]][s,], Obs = cond_y[[s]][[i]][[ind]], n_cat = q_emiss[ind], mu_int_bar1 = c(t(emiss_mu_int_bar[[i]][[q]]) %*% xx[[1 + ind]][s,]), V_int1 = emiss_V_int[[i]][[q]], scalar = emiss_scalar[[q]], candcov1 = emiss_candcov_comb)
            emiss[[s]][[ind]][i,]		   <- PD_subj[[s]]$cat_emiss[iter, (sum(start[1:q]) + 1 + (i - 1) * q_emiss[ind]):(sum(start[1:q]) + (i - 1) * q_emiss[ind] + q_emiss[ind])] <- emiss_RWout$prob + .0001
            emiss_naccept[[q]][s, i]	 <- emiss_naccept[[q]][s, i] + emiss_RWout$accept
            emiss_c_int[[i]][[q]][s,] <- emiss_RWout$draw_int
            emiss_int_subj[[s]][[q]][iter, (1 + (i - 1) * (q_emiss[ind] - 1)) : ((q_emiss[ind] - 1) + (i - 1) * (q_emiss[ind] - 1))]	<- emiss_c_int[[i]][[q]][s,]
          }
        }
      }
      if(n_cont > 0){
        # continuous dependent variables:

        for(q in 1:n_cont){
          ind <- which_cont[q]

          # population parameters
          # mean (and regression parameters if covariates) of the Normal emission distribution
          # and it's variance (so the variance between the subject specific means)
          emiss_mu0_n                    <- solve(t(xx[[1 + ind]]) %*% xx[[1 + ind]] + emiss_cont_K0[[q]]) %*% (t(xx[[1 + ind]]) %*% emiss_c_mu[[i]][[q]] + emiss_cont_K0[[q]] %*% emiss_cont_mu0[[q]][,i])
          emiss_a_mu_n                   <- (emiss_cont_K0[[q]] + n_subj) / 2
          emiss_b_mu_n                   <- (emiss_cont_nu[[q]] * emiss_cont_V[[q]][i]) / 2 + (t(emiss_c_mu[[i]][[q]]) %*% emiss_c_mu[[i]][[q]] +
                                                                                       t(emiss_cont_mu0[[q]][,i]) %*% emiss_cont_K0[[q]] %*% emiss_cont_mu0[[q]][,i] -
                                                                                       t(emiss_mu0_n) %*% (t(xx[[1 + ind]]) %*% xx[[1 + ind]] + emiss_cont_K0[[q]]) %*% emiss_mu0_n) / 2
          emiss_V_mu[[i]][[q]]       <- solve(rgamma(1, shape = emiss_a_mu_n, rate = emiss_b_mu_n))
          if(all(dim(emiss_V_mu[[i]][[q]]) == c(1,1))){ # CHECH THIS
            # emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(as.numeric(emiss_V_mu[[i]][[q]]) * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]))) # This step may produce negative values which is not acceptable
            emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + ind] - 1, mean = 0, sd = sqrt(diag(as.numeric(emiss_V_mu[[i]][[q]]) * solve(t(xx[[1 + ind]]) %*% xx[[1 + ind]] + emiss_cont_K0[[q]]))))
          } else {
            # emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(emiss_V_mu[[i]][[q]] * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]])))
            emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + ind] - 1, mean = 0, sd = sqrt(diag(emiss_V_mu[[i]][[q]] * solve(t(xx[[1 + ind]]) %*% xx[[1 + ind]] + emiss_cont_K0[[q]]))))
          }
          # emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(emiss_V_mu[[i]][[q]] * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]])))
          # if(i > 1){
          #   if(emiss_c_mu_bar[[i]][[q]] < emiss_c_mu_bar[[i-1]][[q]]){
          #     label_switch[i,q] <- label_switch[i,q] + 1
          #
          #   }
          # }

          # fixed variance of the continuous dependent emission distributions
          for (s in 1:n_subj){
            ss_subj[s] <- t(matrix(cond_y[[s]][[i]][[ind]] - emiss_c_mu[[i]][[q]][s,1], nrow = 1) %*%
                              matrix(cond_y[[s]][[i]][[ind]] - emiss_c_mu[[i]][[q]][s,1], ncol = 1))
            n_cond_y[s]       <- length(cond_y[[s]][[i]][[ind]])
          }
          emiss_a_resvar_n <- sum(n_cond_y) / 2 + emiss_cont_a0[[q]][i]
          emiss_b_resvar_n <- (sum(ss_subj) + 2 * emiss_cont_b0[[q]][i]) / 2
          emiss_c_V[[i]][[q]] <- emiss_var_bar[[q]][iter, i] <- solve(rgamma(1, shape = emiss_a_resvar_n, rate = emiss_b_resvar_n))

          ### sampling subject specific means for the emission distributions, assuming known mean and var, see Lynch p. 244
          emiss_c_V_subj    <- (emiss_V_mu[[i]][[q]] * emiss_c_V[[i]][[q]]) / (2 * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
          for (s in 1:n_subj){
            emiss_mu0_subj_n  <- (emiss_V_mu[[i]][[q]] * sum(cond_y[[s]][[i]][[ind]]) +  emiss_c_V[[i]][[q]] * c(t(emiss_c_mu_bar[[i]][[q]]) %*% xx[[1 + ind]][s,])) /
              (n_cond_y[s] * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
            emiss[[s]][[ind]][i,1] <- PD_subj[[s]]$cont_emiss[iter, ((q - 1) * m + i)] <- emiss_c_mu[[i]][[q]][s,1] <- rnorm(1, emiss_mu0_subj_n, sqrt(emiss_c_V_subj))
            emiss[[s]][[ind]][i,2] <- PD_subj[[s]]$cont_emiss[iter, (n_cont * m + (q - 1) * m + i)] <- emiss_c_V[[i]][[q]]
          }
        }
      }
    }

    # End of MCMC iteration, save output values --------
    gamma_int_bar[iter, ]				   	<- unlist(lapply(gamma_mu_int_bar, "[",1,))
    if(nx[1] > 1){
      gamma_cov_bar[iter, ]      	<- unlist(lapply(gamma_mu_int_bar, "[",-1,))
    }
    gamma_prob_bar[iter,]			<- unlist(gamma_mu_prob_bar)

    if(n_cont > 0){
      for(q in 1:n_cont){
        ind <- which_cont[q]
        emiss_mu_bar[[q]][iter, ]	<- as.vector(unlist(lapply(
          lapply(emiss_c_mu_bar, "[[", q), "[",1,)
        ))
        if(nx[1 + ind] > 1){
          emiss_cont_cov_bar[[q]][iter, ]  <- as.vector(unlist(lapply(
            # lapply(emiss_mu_int_bar, "[[", q), "[",-1,)
            lapply(emiss_c_mu_bar, "[[", q), "[",-1,) # CHECK THIS
          ))
        }
        emiss_varmu_bar[[q]][iter,]	<- as.vector(unlist(sapply(emiss_V_mu, "[[", q)))
      }
    }
    if(n_cat > 0){
      for(q in 1:n_cat){
        ind <- which_cat[q]
        emiss_int_bar[[q]][iter, ]	<- as.vector(unlist(lapply(
          lapply(emiss_mu_int_bar, "[[", q), "[",1,)
        ))
        if(nx[1 + ind] > 1){
          emiss_cat_cov_bar[[q]][iter, ]  <- as.vector(unlist(lapply(
            lapply(emiss_mu_int_bar, "[[", q), "[",-1,)
          ))
        }
        emiss_prob_bar[[q]][iter,]	<- as.vector(unlist(sapply(emiss_mu_prob_bar, "[[", q)))
      }
    }

    if(show_progress == TRUE){
      utils::setTxtProgressBar(pb, iter)
    }
  }
  if(show_progress == TRUE){
    close(pb)
  }
  label_switch <- round(label_switch / J * 100, 2)

  # End of function, return output values --------
  ctime = proc.time()[3]
  message(paste("Total time elapsed (hh:mm:ss):", hms(ctime-itime)))
  if(return_path == TRUE){
    if(n_cat * n_cont > 0){
      out <- list(input = list(m = m, n_dep = n_dep , q_emiss = q_emiss, J = J, burn_in = burn_in, data_distr = data_distr,
                               n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj,
                  gamma_int_subj = gamma_int_subj, gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar,
                  gamma_prob_bar = gamma_prob_bar, gamma_naccept = gamma_naccept,
                  emiss_cont_cov_bar = emiss_cont_cov_bar, emiss_cat_cov_bar = emiss_cat_cov_bar, emiss_mu_bar = emiss_mu_bar,
                  emiss_int_subj = emiss_int_subj, emiss_int_bar = emiss_int_bar,
                  emiss_prob_bar = emiss_prob_bar, emiss_naccept = emiss_naccept,
                  emiss_varmu_bar = emiss_varmu_bar, emiss_var_bar = emiss_var_bar,
                  sample_path = sample_path, label_switch = label_switch)
    } else if (n_cat == n_dep){
      out <- list(input = list(m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,burn_in = burn_in,  data_distr = data_distr,
                               n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj,
                  gamma_int_subj = gamma_int_subj, gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar,
                  gamma_prob_bar = gamma_prob_bar, gamma_naccept = gamma_naccept,
                  emiss_int_subj = emiss_int_subj, emiss_int_bar = emiss_int_bar,
                  emiss_cat_cov_bar = emiss_cat_cov_bar, emiss_prob_bar = emiss_prob_bar,
                  emiss_naccept = emiss_naccept,
                  sample_path = sample_path)
    } else if (n_cont == n_dep){
      out <- list(input = list(m = m, n_dep = n_dep, J = J,
                               burn_in = burn_in, , data_distr = data_distr, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj, gamma_int_subj = gamma_int_subj,
                  gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar,
                  emiss_cont_cov_bar = emiss_cont_cov_bar, gamma_prob_bar = gamma_prob_bar,
                  emiss_mu_bar = emiss_mu_bar, gamma_naccept = gamma_naccept,
                  emiss_varmu_bar = emiss_varmu_bar, emiss_var_bar = emiss_var_bar,
                  sample_path = sample_path, label_switch = label_switch)
    }
  } else {
    if(n_cat * n_cont > 0){
      out <- list(input = list(m = m, n_dep = n_dep , q_emiss = q_emiss, J = J, burn_in = burn_in, data_distr = data_distr,
                               n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj,
                  gamma_int_subj = gamma_int_subj, gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar,
                  gamma_prob_bar = gamma_prob_bar, gamma_naccept = gamma_naccept,
                  emiss_cont_cov_bar = emiss_cont_cov_bar, emiss_cat_cov_bar = emiss_cat_cov_bar, emiss_mu_bar = emiss_mu_bar,
                  emiss_int_subj = emiss_int_subj, emiss_int_bar = emiss_int_bar,
                  emiss_prob_bar = emiss_prob_bar, emiss_naccept = emiss_naccept,
                  emiss_varmu_bar = emiss_varmu_bar, emiss_var_bar = emiss_var_bar,
                  label_switch = label_switch)
    } else if (n_cat == n_dep){
      out <- list(input = list(m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,burn_in = burn_in,  data_distr = data_distr,
                               n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj,
                  gamma_int_subj = gamma_int_subj, gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar,
                  gamma_prob_bar = gamma_prob_bar, gamma_naccept = gamma_naccept,
                  emiss_int_subj = emiss_int_subj, emiss_int_bar = emiss_int_bar,
                  emiss_cat_cov_bar = emiss_cat_cov_bar, emiss_prob_bar = emiss_prob_bar,
                  emiss_naccept = emiss_naccept)
    } else if (n_cont == n_dep){
      out <- list(input = list(m = m, n_dep = n_dep, J = J,
                               burn_in = burn_in, data_distr = data_distr, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                  PD_subj = PD_subj, gamma_int_subj = gamma_int_subj,
                  gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar,
                  emiss_cont_cov_bar = emiss_cont_cov_bar, gamma_prob_bar = gamma_prob_bar,
                  emiss_mu_bar = emiss_mu_bar, gamma_naccept = gamma_naccept,
                  emiss_varmu_bar = emiss_varmu_bar, emiss_var_bar = emiss_var_bar,
                  label_switch = label_switch)
    }
  }
  class(out) <- append(class(out), "mHMM_vary")
  return(out)
}
