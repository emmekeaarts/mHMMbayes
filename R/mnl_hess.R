#' @keywords internal
# Obtain mnl -Expected[Hessian]  for intercept only model, bassed on P.Rossi 2004
mnlHess_int <- function(int, Obs, n_cat){
  n_Obs 	<- length(Obs)
  betas   <- c(0, int)
  prob    <- exp(betas) / sum(exp(betas))
  Hess    <- (diag(prob[-1]) - outer(prob[-1],prob[-1])) * n_Obs
  return(Hess)
}
