#' @keywords internal
"_PACKAGE"

# Package imports
# not sure if all functions given below for packages are actually still used, check!
## usethis namespace: start
#' @importFrom mvtnorm dmvnorm rmvnorm dmvt rmvt
#' @importFrom MCMCpack rdirichlet rwish
#' @importFrom stats optim rnorm runif median na.omit quantile
# to include bib references in help file
#' @importFrom Rdpack reprompt
# for RCpp
#' @useDynLib mHMMbayes, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
