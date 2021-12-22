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
#' @seealso \code{\link{mHMM}}, \code{\link{mHMM_cont}} and
#'   \code{\link{mHMM_vary}} for fitting the multilevel hidden Markov model.
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
  if (!is.mHMM(object) & !is.mHMM_cont(object) & !is.mHMM_vary(object)){
    stop("The input object used should either be from the class mHMM, mHMM_cont or mHMM_vary obtained by using their respective functions.")
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
  if(is.mHMM(object)){
    q_emiss <- input$q_emiss
  }
  if(is.mHMM_vary(object)){
    data_distr <- input$data_distr
    n_cat      <- sum(data_distr %in% 'categorical')
    n_cont     <- sum(data_distr %in% 'continuous')
    which_cat  <- which(data_distr %in% 'categorical')
    which_cont <- which(data_distr %in% 'continuous')
    if(n_cat > 0){
      q_emiss <- input$q_emiss
    }
  }
  n_dep   <- input$n_dep
  if (level == "group"){
    if(is.mHMM(object)){
      est <- lapply(q_emiss, dif_matrix, rows = m)
      for(j in 1:n_dep){
        colnames(est[[j]]) <- paste("Category", 1:q_emiss[j])
        rownames(est[[j]]) <- paste("State", 1:m)
      }
      names(est) <- dep_labels
      for(j in 1: n_dep){
        est[[j]][] <- matrix(round(apply(object$emiss_prob_bar[[j]][((burn_in + 1): J),], 2, median),3),
                             byrow = TRUE, ncol = q_emiss[j], nrow = m)
      }
    } else if (is.mHMM_vary(object)){
      est <- lapply(q_emiss, dif_matrix, rows = m)
      for(j in which_cat){
        colnames(est[[j]]) <- paste("Category", 1:q_emiss[j])
        rownames(est[[j]]) <- paste("State", 1:m)
      }
      names(est) <- dep_labels
      if(n_cat > 0){
        for(q in 1:n_cat){
          ind <- which_cat[q]
          est[[ind]][] <- matrix(round(apply(object$emiss_prob_bar[[q]][((burn_in + 1): J),], 2, median),3),
                               byrow = TRUE, ncol = q_emiss[ind], nrow = m)
        }
      }
      if(n_cont > 0){
        for(q in 1:n_cont){
          ind <- which_cont[q]
          est[[ind]] <- matrix(round(c(apply(object$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, median),
                                       apply(object$emiss_var_bar[[q]][((burn_in + 1): J),], 2, median)),3),
                               ncol = 2, nrow = m, dimnames = list(paste("State", 1:m), c("Mean", "Variance")))
        }
      }
    } else if (is.mHMM_cont(object)){
      est <- rep(list(matrix(, nrow = m, ncol = 2, dimnames = list(paste("State", 1:m), c("Mean", "Variance")))), n_dep)
      names(est) <- dep_labels
      for(j in 1:n_dep){
        est[[j]][] <-  matrix(round(c(apply(object$emiss_mu_bar[[j]][((burn_in + 1): J),], 2, median), apply(object$emiss_var_bar[[j]][((burn_in + 1): J),], 2, median)),3), ncol = 2, nrow = m)
      }
    }
    est_emiss <- est
  }
  if (level == "subject"){
    est_emiss <- vector("list", n_dep)
    names(est_emiss) <- dep_labels
    if(is.mHMM(object)){
      for(j in 1:n_dep){
        est <- matrix(,ncol = q_emiss[j], nrow = m,
                      dimnames = list( paste("State", 1:m),
                                       paste("Category", 1:q_emiss[j])))
        est_emiss[[j]] <- rep(list(est), n_subj)
        names(est_emiss[[j]]) <- paste("Subject", 1:n_subj)
      }
      start <- c(0, q_emiss * m)
      for(i in 1:n_subj){
        for(j in 1:n_dep){
          est_emiss[[j]][[i]][] <- matrix(round(apply(object$PD_subj[[i]][burn_in:J, (sum(start[1:j]) + 1) : sum(start[1:(j+1)])], 2, median), 3),
                                          byrow = TRUE, ncol = q_emiss[j], nrow = m)
        }
      }
    } else if (is.mHMM_vary(object)){
      if(n_cat > 0){
        start <- c(0, q_emiss * m)
        for(q in 1:n_cat){
          ind <- which_cat[q]
          est <- matrix(,ncol = q_emiss[ind], nrow = m,
                        dimnames = list( paste("State", 1:m),
                                         paste("Category", 1:q_emiss[ind])))
          est_emiss[[ind]] <- rep(list(est), n_subj)
          names(est_emiss[[ind]]) <- paste("Subject", 1:n_subj)
          for(i in 1:n_subj){
            est_emiss[[ind]][[i]][] <- matrix(round(apply(object$PD_subj[[i]]$cat_emiss[burn_in:J, (sum(start[1:q]) + 1) : sum(start[1:(q+1)])], 2, median), 3),
                                              byrow = TRUE, ncol = q_emiss[ind], nrow = m)
          }
        }
      }
      if(n_cont > 0){
        for(q in 1:n_cont){
          ind <- which_cont[q]
          est_emiss[[ind]] <- rep(list(matrix(, nrow = m, ncol = 2,
                                              dimnames = list(paste("State", 1:m), c("Mean", "Variance")))), n_subj)
          names(est_emiss[[ind]]) <- paste("Subject", 1:n_subj)
          for(i in 1:n_subj){
            est_emiss[[ind]][[i]][] <- matrix(round(c(apply(object$PD_subj[[i]]$cont_emiss[((burn_in + 1): J),((q-1) * m + 1) : ((q-1) * m + m)], 2, median),
                                                      apply(object$PD_subj[[i]]$cont_emiss[((burn_in + 1): J),(n_cont * m + (q-1) * m + 1) : (n_cont * m + (q-1) * m + m)], 2, median)),3), ncol = 2, nrow = m)
          }
        }
      }
    } else if (is.mHMM_cont(object)){
      est_emiss <- rep(list(rep(list(matrix(, nrow = m, ncol = 2, dimnames = list(paste("State", 1:m), c("Mean", "Variance")))), n_subj)), n_dep)
      names(est_emiss) <- dep_labels
      for(j in 1:n_dep){
        names(est_emiss[[j]]) <- paste("Subject", 1:n_subj)
        for(i in 1:n_subj){
          est_emiss[[j]][[i]][] <- matrix(round(c(apply(object$PD_subj[[i]][((burn_in + 1): J),((j-1) * m + 1) : ((j-1) * m + m)], 2, median), apply(object$PD_subj[[i]][((burn_in + 1): J),(n_dep * m + (j-1) * m + 1) : (n_dep * m + (j-1) * m + m)], 2, median)),3), ncol = 2, nrow = m)
        }
      }
    }
  }
  class(est_emiss) <- append(class(est_emiss), "mHMM_emiss")

  return(est_emiss)
}




