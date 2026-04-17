#' Obtain the emission distribution probabilities for a fitted multilevel HMM
#'
#' \code{obtain_emiss_f} obtains the emission distribution for an object
#' containing a fitted multilevel hidden Markov model with a fixed over subjects emission distribution.
#'
#' @inheritParams obtain_gamma
#'
#' @return \code{obtain_emiss_f} creates an object of the class
#'   \code{mHMM_emiss_f}. The output is a list of matrices with the fixed over
#'   subjects emission distribution for each dependent variable.
#'
#' @seealso \code{\link{mHMM_f}} for fitting the
#'   multilevel hidden Markov model with a fixed emission distribution.
#'
#'
#' @examples
#' ###### Currently no example
#' @export
#'
obtain_emiss_f <- function(object, burn_in = NULL){
  if (!is.mHMM_f(object)){
    stop("The input object used should be from the class mHMM_f, obtained by using the function mHMM_f.")
  }
  if(sum(objects(object$PD_subj[[1]]) %in% "log_likl") != 1){
    stop("The input object is created using an earlier version of the mHMMbayes package. Please re-run the function mHMM with the current package version, or post-process the object using the earlier version of the package.")
  }
  input   <- object$input
  dep_labels <- input$dep_labels
  n_subj  <- input$n_subj
  if (is.null(burn_in)){
    burn_in <- input$burn_in
  }
  J       <- input$J
  if (burn_in >= (J-1)){
    stop(paste("The specified burn in period should be at least 2 points smaller
               compared to the number of iterations J, J =", J))
  }
  m       <- input$m
  data_distr <- input$data_distr
  n_dep   <- input$n_dep

  est <- rep(list(matrix(NA_real_, nrow = m, ncol = 2, dimnames = list(paste("State", 1:m), c("Mean", "SD")))), n_dep)
  names(est) <- dep_labels
  for(q in 1:n_dep){
    est[[q]][] <-  matrix(round(c(apply(matrix(object$emiss_mu[[q]][((burn_in + 1): J),], ncol = m), 2, median), apply(matrix(object$emiss_sd[[q]][((burn_in + 1): J),], ncol = m), 2, median)),3), ncol = 2, nrow = m)
  }
  est_emiss <- est
  class(est_emiss) <- append("mHMM_emiss_f", class(est_emiss))
  return(est_emiss)
}




