#' @export
#'
print.HMM <- function(x, ...){
  input   <- x$input
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  q_emiss <- input$q_emiss
  n_dep   <- input$n_dep
  cat(J, "iterations used in the MCMC algorithm with a burn in of", burn_in, "\n")
  LL      <- median(x$PD[((burn_in + 1): J), (sum(q_emiss * m) + m*m + 1)])
  AIC<-2*(sum((q_emiss-1)*m)+(m-1)*m) - (2*LL)
  cat("\n")
  cat("Median Log likelihood over the iterations:", LL, "\n")
  cat("Medain AIC over the iterations:", AIC, "\n")
  cat("\n")
  cat("Number of states used:", m, "\n")
  cat("\n")
  cat("Number of dependent variables used:", n_dep, "\n")
  cat("\n")
}
