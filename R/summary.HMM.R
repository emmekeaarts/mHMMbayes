#' @export
#'
summary.HMM <- function(object, ...){
  input   <- object$input
  dep_labels <- input$dep_labels
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  q_emiss <- input$q_emiss
  n_dep   <- input$n_dep
  gamma_est <- matrix(round(apply(object$PD[((burn_in + 1): J), (sum(q_emiss * m) + 1): (sum(q_emiss * m) + m*m)], 2, median),3), byrow = TRUE, ncol = m, nrow = m)
  colnames(gamma_est) <- paste("To state", 1:m)
  rownames(gamma_est) <- paste("From state", 1:m)
  cat("State transition probability matrix:", "\n", "\n")
  print(gamma_est)
  cat("\n", "\n")
  cat("Emission distribution for each of the dependent variables:", "\n", "\n")
  EM_est <- vector("list", n_dep)
  names(EM_est) <- dep_labels
  start <- c(0, q_emiss * m)
  for(i in 1:n_dep){
    EM_est[[i]] <- matrix(round(apply(object$PD[((burn_in + 1): J),(sum(start[1:i]) + 1) : sum(start[1:(i+1)])], 2, median),3), byrow = TRUE, ncol = q_emiss[i], nrow = m)
    colnames(EM_est[[i]]) <- paste("Category", 1:q_emiss[i])
    rownames(EM_est[[i]]) <- paste("State", 1:m)
  }
  print(EM_est)
  cat("\n")
}
