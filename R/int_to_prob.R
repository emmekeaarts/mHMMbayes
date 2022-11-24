#' Transforming a set of Multinomial logit regression intercepts to probabilities
#'
#' \code{int_to_prob} transforms a set of Multinomial logit regression
#' intercepts to the corresponding state transition or categorical emission
#' observation probabilities. Note that the first state or category is assumed
#' to be the reference category, hence no intercept is to specified for the
#' first state or category.
#'
#' Designed to ease the specification of informative hyper-prior values for the
#' mean intercepts of the transition probability matrix gamma and categorical
#' emission distribution(s) of the multilevel hidden Markov model through the
#' functions \code{\link{prior_gamma}} and \code{\link{prior_emiss_cat}}. No
#' check is performed on correct specifications of the dimensions.
#'
#' @param int_matrix A matrix with (number of states OR categories - 1) columns
#'   and number of rows to be determined by the user. For obtaining the set of
#'   probabilities of the complete transition probability matrix gamma or
#'   categorical emission distribution matrix, the number of rows equals the
#'   number of states \code{m}. The first state / category is assumed to be the
#'   reference category, no intercept is to be specified for this first
#'   category.
#'
#' @return \code{int_to_prob} returns a matrix containing probabilities with
#'   each row summing to one, with the number of columns equal to the number of
#'   states / categories and the number of rows equal to the number rows
#'   specified in the input matrix.
#'
#' @seealso \code{\link{prob_to_int}} for transforming a set of probabilities to
#'   a set of Multinomial logit regression intercepts, \code{\link{prior_gamma}}
#'   and \code{\link{prior_emiss_cat}} for specifying informative hyper-priors
#'   for the the multilevel hidden Markov model and \code{\link{mHMM}} to fit a
#'   multilevel hidden Markov model.
#'
#' @examples
#'
#' # example for transition probability matrix gamma with 3 states
#' m <- 3
#' gamma_int <- matrix(c(-1, -1,
#'                        3,  0,
#'                        0,  2), ncol = m-1, nrow = m, byrow = TRUE)
#' gamma_prob <- int_to_prob(gamma_int)
#' gamma_prob
#'
#' @export

int_to_prob <- function(int_matrix) {
  if(!is.matrix(int_matrix)){
    stop("int_matrix should be a matrix")
  }
  prob_matrix <- matrix(nrow = nrow(int_matrix), ncol = ncol(int_matrix) + 1)
  for(r in 1:nrow(int_matrix)){
    exp_int_matrix 	<- matrix(exp(c(0, int_matrix[r,])), nrow  = 1)
    prob_matrix[r,] <- exp_int_matrix / as.vector(exp_int_matrix %*% c(rep(1, (dim(exp_int_matrix)[2]))))
    }
  return(round(prob_matrix,4))
}

#' Transforming a set of probabilities to Multinomial logit regression intercepts
#'
#' \code{prob_to_int} transforms a set of state transition or categorical
#' emission observation probabilities to the corresponding Multinomial logit
#' regression intercepts. Note that the first category is assumed to be the
#' reference category, hence no intercept is returned for the first state or
#' category.
#'
#' Designed to ease the specification of informative hyper-prior values for the
#' mean intercepts of the transition probability matrix gamma and categorical
#' emission distribution(s) of the multilevel hidden Markov model through the
#' functions \code{\link{prior_gamma}} and \code{\link{prior_emiss_cat}}. No
#' check is performed on correct specifications of the dimensions.
#'
#' @param prob_matrix A matrix with number of states OR categories columns and
#'   number of rows to be determined by the user, with rows summing to one. For
#'   obtaining the set of Multinomial logit regression intercepts of the
#'   complete transition probability matrix gamma or categorical emission
#'   distribution matrix, the number of rows equals the number of states
#'   \code{m}.
#'
#' @return \code{prob_to_int} returns a matrix containing Multinomial logit
#'   regression intercepts, with the number of columns equal to (number of
#'   states or categories - 1) and the number of rows equal to the number rows
#'   specified in the input matrix. The first state / category is assumed to be
#'   the reference category, hence no intercept is returned for this first
#'   category.
#'
#' @seealso \code{\link{int_to_prob}} for transforming a set of Multinomial
#'   logit regression intercepts to a probabilities, \code{\link{prior_gamma}}
#'   and \code{\link{prior_emiss_cat}} for specifying informative hyper-priors
#'   for the the multilevel hidden Markov model and \code{\link{mHMM}} to fit a
#'   multilevel hidden Markov model.
#'
#' @examples
#'
#' # example for transition probability matrix gamma with 3 states
#' m <- 3
#' gamma_prob <- matrix(c(0.6, 0.2, 0.2,
#'                        0.1, 0.8, 0.1,
#'                        0.1, 0.1, 0.8), ncol = m, nrow = m, byrow = TRUE)
#' gamma_int <- prob_to_int(gamma_prob)
#' gamma_int
#'
#' @export

prob_to_int <- function(prob_matrix){
  if(!is.matrix(prob_matrix)){
    stop("prob_matrix should be a matrix")
  }
  prob_matrix <- prob_matrix + 0.00001
  int_matrix <- matrix(NA, nrow(prob_matrix), ncol(prob_matrix)-1)
  sum_exp <- numeric(nrow(prob_matrix))
  for(r in 1:nrow(prob_matrix)){
    sum_exp[r] <- (1/prob_matrix[r,1]) - 1
    for(cr in 2:ncol(prob_matrix)){
      #for every int_matrix except the first collumn (e.g. int_matrix12 <- log(y12/y11-y12))
      int_matrix[r,(cr-1)] <- log(prob_matrix[r,cr]*(1+sum_exp[r]))
    }
  }
  return(round(int_matrix,4))
}
