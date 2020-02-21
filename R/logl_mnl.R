#' @keywords internal
# Evaluate loglikelihood for intercept only MNL, bassed on P. Rossi 2004
llmnl_int <- function(int, Obs, n_cat) {
  n_Obs <- length(Obs)
  betas <- c(0, int)
  return(sum(betas[Obs]) - log(sum(exp(betas))) * n_Obs)
}

#' @keywords internal
# Obtain fractional log likelihood for multinomial intercept only model, bassed on P. Rossi 2004
llmnl_int_frac <- function(int, Obs, n_cat, pooled_likel, w, wgt){
  return((1 - w) * llmnl_int(int = int, Obs = Obs, n_cat = n_cat) + w * wgt * pooled_likel)
}
