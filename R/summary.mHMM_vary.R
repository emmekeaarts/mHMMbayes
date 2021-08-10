#' @export
#'
summary.mHMM_vary <- function(object, ...){
  input   <- object$input
  dep_labels <- input$dep_labels
  n_subj  <- input$n_subj
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  n_dep   <- input$n_dep
  data_distr <- input$data_distr
  n_cat      <- sum(data_distr %in% 'categorical')
  n_cont     <- sum(data_distr %in% 'continuous')
  which_cat  <- which(data_distr %in% 'categorical')
  which_cont <- which(data_distr %in% 'continuous')
  if(n_cat > 0){
    q_emiss <- input$q_emiss
  }

  gamma_pop <- matrix(round(apply(object$gamma_prob_bar[((burn_in + 1): J),], 2, median),3), byrow = TRUE, ncol = m, nrow = m)
  colnames(gamma_pop) <- paste("To state", 1:m)
  rownames(gamma_pop) <- paste("From state", 1:m)
  cat("State transition probability matrix","\n",  "(at the group level):", "\n", "\n")
  print(gamma_pop)
  cat("\n", "\n")
  cat("Emission distribution for each of the dependent variables","\n",  "(at the group level):", "\n", "\n")
  EM_pop <- vector("list", n_dep)
  names(EM_pop) <- dep_labels

  if(n_cat > 0){
    for(q in 1:n_cat){
      ind <- which_cat[q]
      EM_pop[[ind]] <- matrix(round(apply(object$emiss_prob_bar[[q]][((burn_in + 1): J),], 2, median),3), byrow = TRUE, ncol = q_emiss[ind], nrow = m)
      colnames(EM_pop[[ind]]) <- paste("Category", 1:q_emiss[ind])
      rownames(EM_pop[[ind]]) <- paste("State", 1:m)
    }
  }
  if(n_cont > 0){
    for(q in 1:n_cont){
      ind <- which_cont[q]
      EM_pop[[ind]] <- matrix(round(c(apply(object$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, median), apply(object$emiss_var_bar[[q]][((burn_in + 1): J),], 2, median)),3), ncol = 2, nrow = m)
      colnames(EM_pop[[ind]]) <- c("Mean", "Variance")
      rownames(EM_pop[[ind]]) <- paste("State", 1:m)
    }
  }
  print(EM_pop)
  cat("\n")
}
