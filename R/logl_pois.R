#' @keywords internal
# Obtain poisson log likelihood
llpois <- function(lambda, Obs){
  ll_pois <- sum(stats::dpois(Obs, lambda, log = TRUE))
  return(ll_pois)
}

#' @keywords internal
# Obtain fractional poisson log likelihood
llpois_frac_log <- function(par, Obs, n_cat, pooled_likel, w, wgt){
  return((1 - w) * llpois(lambda = exp(par), Obs = Obs) + w * wgt * pooled_likel)
}
