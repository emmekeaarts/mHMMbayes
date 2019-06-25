#' mHMMbayes: A package for multilevel hidden Markov models using Bayesian estimation.
#'
#'   With the \code{R} package \code{mHMMbayes} you can fit multilevel hidden
#'   Markov models. The multilevel hidden Markov model (HMM) is a generalization
#'   of the well-known hidden Markov model, tailored to accomodate (intense)
#'   longitudinal data of multiple individuals simultaneously. Using a
#'   multilevel framework, we allow for heterogeneity in the model parameters
#'   (transition probability matrix and conditional distribution), while
#'   estimating one overall HMM. The model has a great potential of application
#'   in  many fields, such as the social sciences and medicine. The model can be
#'   fitted on multivariate data with a catagorical  distribution, and include
#'   individual level covariates (allowing for e.g., group comparisons on model
#'   parameters). Parameters are estimated using Bayesian estimation utilizing
#'   the forward-backward recursion within a hybrid Metropolis within Gibbs
#'   sampler.
#'
#'   The \code{mHMMbayes} package provides three main functions: \code{\link{mHMM_mnl}}
#'   , \code{\link{sim_mHMM}} and \code{\link{vit_mHMM}}, described below. For a more
#'   elaborate guide to the package \code{mHMMbayes}, see the tutorial-mhmm
#'   vignette: \code{vignette("tutorial-mhmm", package = "mHMMbayes")} . For extensive
#'   information on the estimation of the parameters in the package, see the estimation-mhmm
#'   vignette: \code{vignette("estimation-mhmm", package = "mHMMbayes")}.
#'
#' @section \code{mHMM_mnl}:
#' The function \code{mHMM_mnl} is used to analyses (intense longitudinal) data
#' from multiple subjects using a multilevel hidden Markov model. By using a
#' multilevel framework, one general 'population' HMM is estimated, while
#' heterogeneity between subjects is accommodated. The function can handle
#' covariates at the subject level (unlimited number), uses a hybrid metropolis
#' within gibs sampler, and performs the forward backward algorithm for all
#' subjects in a sequential manner. Can handle varying observation length over
#' subjects.
#'
#' @section \code{sim_mHMM}:
#' The function \code{sim_mHMM} is used to simulate data for multiple subjects,
#' for which the data have categorical observations that follow a hidden Markov
#' model (HMM) with an multilevel structure. The multilevel structure implies
#' that each subject is allowed to have it's own set of parameters, and that the
#' parameters at the subject level (level 1) are tied together by a population
#' distribution at level 2 for each of the corresponding parameters. The shape
#' of the population distribution for each of the parameters is a normal (i.e.,
#' Gaussian) distribution. In addition to (natural and/or unexplained)
#' heterogeneity between subjects, the subjects parameters can also depend on a
#' (set of) covariate(s).
#'
#' @section \code{vit_mHMM}:
#' The function \code{vit_mHMM} is used to obtain the most likely hidden state
#' sequence for each subject, given the data and the subject specific parameter
#' estimates. The function does this by utilizing the Viterbi algorithm.
#'
# Package imports
# not sure if all functions given below for packages are actually still used, check!
#' @importFrom mvtnorm dmvnorm rmvnorm dmvt rmvt
#' @importFrom MCMCpack rdirichlet rwish
#' @importFrom stats optim rnorm runif median
#' @importFrom alluvial alluvial
#'
#'
#' @docType package
#' @name mHMMbayes
NULL
