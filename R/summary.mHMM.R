#' @export
#'
summary.mHMM <- function(object, ...){
  if(sum(objects(object$PD_subj[[1]]) %in% "log_likl") != 1){
    stop("The input object is created using an earlier version of the mHMMbayes package. Please re-run the function mHMM with the current package version, or post-process the object using the earlier version of the package.")
  }
  input   <- object$input
  dep_labels <- input$dep_labels
  n_subj  <- input$n_subj
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  data_distr <- input$data_distr
  if(data_distr == 'categorical'){
    q_emiss <- input$q_emiss
  }
  n_dep   <- input$n_dep
  if(m > 1){
    gamma_int <- matrix(apply(object$gamma_int_bar[((burn_in + 1): J),], 2, median), byrow = TRUE, ncol = m-1, nrow = m)
    gamma_pop <- round(int_to_prob(gamma_int),3)
  } else if (m == 1){
    gamma_pop <- matrix(1)
  }
  colnames(gamma_pop) <- paste("To state", 1:m)
  rownames(gamma_pop) <- paste("From state", 1:m)
  cat("State transition probability matrix","\n",  "(at the group level):", "\n", "\n")
  print(gamma_pop)
  cat("\n", "\n")
  cat("Emission distribution (", data_distr, ") for each of the dependent variables","\n",  "(at the group level):", "\n", "\n")
  EM_pop <- vector("list", n_dep)
  names(EM_pop) <- dep_labels
  if(data_distr == 'categorical'){
    EM_int <- vector("list", n_dep)
    for(i in 1:n_dep){
      EM_int[[i]] <- matrix(apply(object$emiss_int_bar[[i]][((burn_in + 1): J),], 2, median), byrow = TRUE, ncol = q_emiss[i]-1, nrow = m)
      EM_pop[[i]] <- round(int_to_prob(EM_int[[i]]),3)
      colnames(EM_pop[[i]]) <- paste("Category", 1:q_emiss[i])
      rownames(EM_pop[[i]]) <- paste("State", 1:m)
    }
  } else if (data_distr == 'continuous'){
    for(i in 1:n_dep){
      EM_pop[[i]] <- matrix(round(c(apply(matrix(object$emiss_mu_bar[[i]][((burn_in + 1): J),], ncol = m), 2, median),
                                    apply(matrix(object$emiss_sd_bar[[i]][((burn_in + 1): J),], ncol = m), 2, median)),
                                  3), ncol = 2, nrow = m)
      colnames(EM_pop[[i]]) <- c("Mean", "SD")
      rownames(EM_pop[[i]]) <- paste("State", 1:m)
    }
  } else if (data_distr == 'count'){
    for(i in 1:n_dep){
      EM_pop[[i]] <- matrix(round(
        apply(matrix(object$emiss_mu_bar[[i]][((burn_in + 1): J),], ncol = m), 2, median),
        3), ncol = 1, nrow = m)
      colnames(EM_pop[[i]]) <- c("Mean")
      rownames(EM_pop[[i]]) <- paste("State", 1:m)
    }
  }
  print(EM_pop)
  cat("\n")
}
