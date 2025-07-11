#' Obtain the transition probabilities gamma for a fitted multilevel HMM
#'
#' \code{obtain_gamma} obtains the transition probability matrix (TPM) for an
#'   object containing a fitted multilevel hidden Markov model. The TPM can be
#'   obtained either at the group level, i.e., representing the average transition
#'   probability matrix over all subjects, or at the subject level, returning the
#'   transition probability matrices for each subject.
#'
#' @param object An object of class \code{mHMM}, generated by the function
#'   \code{\link{mHMM}}.
#' @param level String specifying if the returned transition probability matrix
#'   gamma should be at the group level (\code{level = "group"}), i.e.,
#'   representing the average transition probability matrix over all subjects,
#'   or at the subject level (\code{level = "subject"}).
#' @param burn_in An integer which specifies the number of iterations to discard
#'   when obtaining the model parameter summary statistics. When left
#'   unspecified (\code{burn_in = NULL}), the burn in period specified when
#'   creating the \code{mHMM} object will be used.
#'
#' @return \code{obtain_gamma} creates an object of the class \code{mHMM_gamma}.
#'   This object can be directly plotted using the function
#'   \code{plot.mHMM_gamma()}, or simply \code{plot()}. Depending on the
#'   specification at the input variable \code{level}, the output is either a
#'   matrix with the transition probabilities at the group level (if \code{level
#'   = "group"}), or a list of matrices (with the number of elements equal to
#'   the number of subjects analyzed, if \code{level = 'subject'}), where each
#'   matrix in the list represents a subject specific transition probability
#'   matrix.
#'
#' @seealso \code{\link{mHMM}} for fitting the multilevel hidden Markov model,
#'   and \code{\link{plot.mHMM_gamma}} for plotting the obtained transition
#'   probabilities.
#'
#'
#'
#' @examples
#' ###### Example on package data
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
#' out_2st <- mHMM(s_data = nonverbal,
#'                 gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
#'                 start_val = c(list(start_TM), start_EM),
#'                 mcmc = list(J = 11, burn_in = 5))
#'
#' out_2st
#' summary(out_2st)
#'
#' # obtaining the transition probabilities at the group and subject level
#' obtain_gamma(out_2st, level = "group")
#' obtain_gamma(out_2st, level = "subject")
#'
#' }
#'
#' @export
#'
obtain_gamma <- function(object, level = "group", burn_in = NULL){
  if (!is.mHMM(object)){
    stop("The input object used should be from the class mHMM, obtained by using the function mHMM.")
  }
  if (level != "group" & level != "subject"){
    stop("The specification at the input variable -level- should be set to either group or subject")
  }
  input   <- object$input
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
  if(m == 1){
    stop("obtain_gamma cannot be used when the number of states m equals 1; it's obsolete (the transition probability matrix always equals 1).")
  }
  n_dep   <- input$n_dep
  est <- matrix(, ncol = m, nrow = m)
  colnames(est) <- paste("To state", 1:m)
  rownames(est) <- paste("From state", 1:m)
  if (level == "group"){
    gamma_int <- matrix(apply(object$gamma_int_bar[((burn_in + 1): J),], 2, median), byrow = TRUE, ncol = m-1, nrow = m)
    est[]     <- round(int_to_prob(gamma_int),3)
    est_gamma <- est
  }
  if (level == "subject"){
  est_gamma <- rep(list(est), n_subj)
  names(est_gamma) <- paste("Subject", 1:n_subj)
   for(i in 1:n_subj){
     gamma_int <- matrix(apply(object$gamma_int_subj[[i]][((burn_in + 1): J),], 2, median), byrow = TRUE, ncol = m-1, nrow = m)
     est_gamma[[i]][] <- round(int_to_prob(gamma_int),3)
   }
  }
  class(est_gamma) <- append(class(est_gamma), "mHMM_gamma")
  return(est_gamma)
}
