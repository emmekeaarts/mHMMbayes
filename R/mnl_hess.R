#' @keywords internal
mnlHess_int <- function(int, Obs, n_cat){
  n_Obs 	<- length(Obs)
  betas   <- matrix(c(0, int), byrow = T, ncol = n_cat)
  prob    <- exp(betas) / sum(exp(betas))
  Hess    <- (diag(x = prob[-1], nrow = n_cat-1) - outer(prob[-1],prob[-1])) * n_Obs
  return(Hess)
}
