#' @export
#'
print.mHMM_vary <- function(x, ...){
  input   <- x$input
  n_subj  <- input$n_subj
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  n_dep   <- input$n_dep
  data_distr <- input$data_distr

  cat("Number of subjects:", n_subj, "\n")
  cat("\n")
  cat(J, "iterations used in the MCMC algorithm with a burn in of", burn_in, "\n")
  LL      <- numeric(n_subj)
  for(i in 1:n_subj){
    LL[i] <- median(x$PD_subj[[i]]$log_likl[((burn_in + 1): J), 1])
  }
  AIC<-2*(m * n_dep * 2 + (m - 1) * m) - (2 * LL)
  cat("Average Log likelihood over all subjects:", mean(LL), "\n")
  cat("Average AIC over all subjects:", mean(AIC), "\n")
  cat("\n")
  cat("Number of states used:", m, "\n")
  cat("\n")
  cat("Number of dependent variables used:", n_dep, "\n")
  cat("\n")
  cat("Distribution(s) of the", n_dep, "dependent variable(s) used:", "\n", data_distr, "\n")
  cat("\n")
}
