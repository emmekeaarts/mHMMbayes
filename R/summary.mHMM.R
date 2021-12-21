#' @export
#'
summary.mHMM <- function(object, ...){
  input   <- object$input
  dep_labels <- input$dep_labels
  n_subj  <- input$n_subj
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  q_emiss <- input$q_emiss
  n_dep   <- input$n_dep
  gamma_int_est <- matrix(apply(object$gamma_int_bar[((burn_in + 1): J),], 2, median), byrow = TRUE, ncol = m-1, nrow = m)
  gamma_pop<-round(int_to_prob(gamma_int_est),3)
  colnames(gamma_pop) <- paste("To state", 1:m)
  rownames(gamma_pop) <- paste("From state", 1:m)
  cat("State transition probability matrix","\n",  "(at the group level):", "\n", "\n")
  print(gamma_pop)
  cat("\n", "\n")
  cat("Emission distribution for each of the dependent variables","\n",  "(at the group level):", "\n", "\n")
  EM_int_est<- vector("list", n_dep)
  EM_pop <- vector("list", n_dep)

  names(EM_pop) <- dep_labels
  for(i in 1:n_dep){
    i=4
    EM_int_est[[i]] <- matrix(apply(object$emiss_int_bar[[i]][((burn_in + 1): J),], 2, median), byrow = TRUE, ncol = q_emiss[i]-1, nrow = m)
    EM_pop[[i]]<-round(int_to_prob(EM_int_est[[i]]),3)

    colnames(EM_pop[[i]]) <- paste("Category", 1:q_emiss[i])
    rownames(EM_pop[[i]]) <- paste("State", 1:m)
  }
  print(EM_pop)
  cat("\n")
}
