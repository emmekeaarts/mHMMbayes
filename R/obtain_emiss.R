#' Obtain the emission distribution probabilities for a fitted multilevel HMM
#'
#' \code{obtain_emiss} obtains the emission distribution for an object
#' containing a fitted multilevel hidden Markov model, either at the group
#' level, i.e., representing the average emission distribution over all
#' subjects, or at the subject level, returning the emission distribution for
#' each subject.
#'
#' @inheritParams obtain_gamma
#'
#' @return \code{obtain_emiss} creates an object of the class \code{mHMM_emiss}.
#'   Depending on the specification at the input variable \code{level}, the
#'   output is either a list of matrices with the emission distribution at the
#'   group level (if \code{level = "group"}) for each dependent variable, or a
#'   list of lists, where for each dependent variable a list is returned with
#'   the number of elements equal to the number of subjects analyzed (if
#'   \code{level = 'subject'}). In the latter scenario, each matrix in the lower
#'   level list represents the subject specific emission distribution for a specific
#'   dependent variable.
#'
#' @seealso \code{\link{mHMM}} for fitting the
#'   multilevel hidden Markov model.
#'
#'
#' @examples
#' ###### Example on package data, see ?nonverbal
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
  if(sum(objects(object$PD_subj[[1]]) %in% "log_likl") != 1){
    stop("The input object is created using an earlier version of the mHMMbayes package. Please re-run the function mHMM with the current package version, or post-process the object using the earlier version of the package.")
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
  data_distr <- input$data_distr
  if(data_distr == 'categorical'){
    q_emiss <- input$q_emiss
  }
  n_dep   <- input$n_dep

  if (level == "group"){
    if(data_distr == 'categorical'){
      est_int <- est <- vector("list", n_dep)
      names(est) <- dep_labels
      for(q in 1:n_dep){
        est_int[[q]] <- matrix(apply(matrix(object$emiss_int_bar[[q]][((burn_in + 1): J),], ncol = (q_emiss[q]-1) * m), 2, median), byrow = TRUE, ncol = q_emiss[q]-1, nrow = m)
        est[[q]] <- round(int_to_prob(est_int[[q]]),3)
        colnames(est[[q]]) <- paste("Category", 1:q_emiss[q])
        rownames(est[[q]]) <- paste("State", 1:m)
      }
    } else if (data_distr == 'continuous'){
      est <- rep(list(matrix(NA_real_, nrow = m, ncol = 2, dimnames = list(paste("State", 1:m), c("Mean", "SD")))), n_dep)
      names(est) <- dep_labels
      for(q in 1:n_dep){
        est[[q]][] <-  matrix(round(c(apply(matrix(object$emiss_mu_bar[[q]][((burn_in + 1): J),], ncol = m), 2, median), apply(matrix(object$emiss_sd_bar[[q]][((burn_in + 1): J),], ncol = m), 2, median)),3), ncol = 2, nrow = m)
      }
    } else if (data_distr == 'count'){
      est <- rep(list(matrix(NA_real_, nrow = m, ncol = 1, dimnames = list(paste("State", 1:m), "Mean"))), n_dep)
      names(est) <- dep_labels
      for(q in 1:n_dep){
        est[[q]][] <-  matrix(round(apply(matrix(object$emiss_mu_bar[[q]][((burn_in + 1): J),], ncol = m), 2, median),3), ncol = 1, nrow = m)
      }
    }

  est_emiss <- est
  }
  if (level == "subject"){
    if(data_distr == 'categorical'){
      est_emiss <- est_emiss_int <- vector("list", n_dep)
      names(est_emiss) <- dep_labels
      for(q in 1:n_dep){
       est <- matrix(NA_real_, ncol = q_emiss[q], nrow = m)
       colnames(est) <- paste("Category", 1:q_emiss[q])
       rownames(est) <- paste("State", 1:m)
       est_emiss[[q]] <- rep(list(est), n_subj)
       names(est_emiss[[q]]) <- paste("Subject", 1:n_subj)
       est_int <-  matrix(NA_real_, ncol = q_emiss[q] - 1, nrow = m)
       est_emiss_int[[q]] <- rep(list(est_int), n_subj)
     }
      for(s in 1:n_subj){
      for(q in 1:n_dep){
          est_emiss_int[[q]][[s]][] <- matrix(apply(matrix(object$emiss_int_subj[[s]][[q]][burn_in:J, ], ncol = (q_emiss[q]-1) * m), 2, median),
                                              byrow = TRUE, ncol = q_emiss[q]-1, nrow = m)
          est_emiss[[q]][[s]][] <- round(int_to_prob(est_emiss_int[[q]][[s]]),3)
        }
      }
    } else if (data_distr == 'continuous'){
      est_emiss <- rep(list(rep(list(matrix(NA_real_, nrow = m, ncol = 2, dimnames = list(paste("State", 1:m), c("Mean", "SD")))), n_subj)), n_dep)
      names(est_emiss) <- dep_labels
      for(q in 1:n_dep){
        names(est_emiss[[q]]) <- paste("Subject", 1:n_subj)
        for(s in 1:n_subj){
          est_emiss[[q]][[s]][] <- matrix(round(c(apply(matrix(object$PD_subj[[s]]$cont_emiss[((burn_in + 1): J),((q-1) * m + 1):(q * m)], ncol = m), 2, median),
                                                  apply(matrix(object$PD_subj[[s]]$cont_emiss[((burn_in + 1): J), (n_dep * m + (q-1) * m + 1):(n_dep * m + q * m)], ncol = m), 2, median)),
                                                3), ncol = 2, nrow = m)
        }
      }
    } else if (data_distr == 'count'){
      est_emiss <- rep(list(rep(list(matrix(NA_real_, nrow = m, ncol = 1, dimnames = list(paste("State", 1:m), c("Mean")))), n_subj)), n_dep)
      names(est_emiss) <- dep_labels
      for(q in 1:n_dep){
        names(est_emiss[[q]]) <- paste("Subject", 1:n_subj)
        for(s in 1:n_subj){
          est_emiss[[q]][[s]][] <- matrix(round(c(apply(matrix(object$PD_subj[[s]]$count_emiss[((burn_in + 1): J),((q-1) * m + 1):(q * m)], ncol = m), 2, median)),3), ncol = 1, nrow = m)
        }
      }
    }
  }
  class(est_emiss) <- append("mHMM_emiss", class(est_emiss))
  return(est_emiss)
}




