#' Obtain the emission distribution probabilities for a fitted multilevel HMM
#'
#' \code{obtain_emiss} obtains the emission distribution for a fitted multilevel
#' hidden Markov model, for either the group level, i.e., representing the
#' average emission distribution over all subjects, or at the subject level,
#' returning the emission distribution for each subject.
#'
#' @inheritParams obtain_gamma
#'
#' @return \code{obtain_emiss} returns the object \code{est_emiss}. Depending on
#'   the specification at the input variable \code{level}, \code{est_emiss} is
#'   either a list of matrices with the emission distribution at the group level
#'   (if \code{level = "group"}) for each dependent variable, or a list of
#'   lists, where for each dependent variable a list is returned with the number
#'   of elements equal to the number of subjects analyzed, if \code{level =
#'   'subject'}). In the latter scenario, each matrix in the lower level list
#'   represents the subject specific emission distribution for a specific
#'   dependent variable.
#'
#' @seealso \code{\link{mHMM}} for fitting the multilevel hidden Markov model.
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
#' # obtaining the emission probabilities at the group and subject level
#' obtain_emiss(out_2st, level = "group")
#' obtain_emiss(out_2st, level = "subject")
#'
#' }
#'
#' @export
#'
obtain_emiss <- function(object, level = "group", burn_in = NULL){
  if (!is.mHMM(object)){
    stop("The input object used should be from the class mHMM, obtained by using the function mHMM.")
  }
  if (level != "group" & level != "subject"){
    stop("The specification at the input variable -level- should be set to either group or subject")
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
  q_emiss <- input$q_emiss
  n_dep   <- input$n_dep
  if (level == "group"){
    est_int <- est <- vector("list", n_dep)
    names(est) <- dep_labels
    for(i in 1:n_dep){
      est_int[[i]] <- matrix(apply(object$emiss_int_bar[[i]][((burn_in + 1): J),], 2, median), byrow = TRUE, ncol = q_emiss[i]-1, nrow = m)
      est[[i]] <- round(int_to_prob(est_int[[i]]),3)
      colnames(est[[i]]) <- paste("Category", 1:q_emiss[i])
      rownames(est[[i]]) <- paste("State", 1:m)
    }

  est_emiss <- est
  }
  if (level == "subject"){
    est_emiss <- est_emiss_int <- vector("list", n_dep)
    names(est_emiss) <- dep_labels
    for(j in 1:n_dep){
      est <- matrix(, ncol = q_emiss[j], nrow = m)
      colnames(est) <- paste("Category", 1:q_emiss[j])
      rownames(est) <- paste("State", 1:m)
      est_emiss[[j]] <- rep(list(est), n_subj)
      names(est_emiss[[j]]) <- paste("Subject", 1:n_subj)
      est_int <-  matrix(, ncol = q_emiss[j] - 1, nrow = m)
      est_emiss_int[[j]] <- rep(list(est_int), n_subj)
    }
    for(i in 1:n_subj){
      for(j in 1:n_dep){
        est_emiss_int[[j]][[i]][] <- matrix(apply(object$emiss_int_subj[[i]][[j]][burn_in:J, ], 2, median),
                                           byrow = TRUE, ncol = q_emiss[j]-1, nrow = m)
        est_emiss[[j]][[i]][] <- round(int_to_prob(est_emiss_int[[j]][[i]]),3)
      }
    }
  }
  return(est_emiss)
}




