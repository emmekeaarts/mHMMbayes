#' @export
#'
summary.mHMM_tr_cont <- function(object, ...){
  input   <- object$input
  dep_labels <- input$dep_labels
  n_subj  <- input$n_subj
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  n_dep   <- input$n_dep
  gamma_int <- matrix(apply(object$gamma_int_bar[((burn_in + 1): J),], 2, median), byrow = TRUE, ncol = m-1, nrow = m)
  gamma_pop <- round(int_to_prob(gamma_int),3)
  colnames(gamma_pop) <- paste("To state", 1:m)
  rownames(gamma_pop) <- paste("From state", 1:m)
  cat("State transition probability matrix","\n",  "(at the group level):", "\n", "\n")
  print(gamma_pop)
  cat("\n", "\n")
  cat("Emission distribution for each of the dependent variables","\n",  "(at the group level):", "\n", "\n")
  EM_pop <- vector("list", n_dep)
  names(EM_pop) <- dep_labels
  for(i in 1:n_dep){
    EM_pop[[i]] <- matrix(round(c(apply(object$emiss_mu_bar[[i]][((burn_in + 1): J),], 2, median), apply(object$emiss_var_bar[[i]][((burn_in + 1): J),], 2, median)),3), ncol = 2, nrow = m)
    colnames(EM_pop[[i]]) <- c("Mean", "Variance")
    rownames(EM_pop[[i]]) <- paste("State", 1:m)
  }
  print(EM_pop)
  cat("\n")
}
