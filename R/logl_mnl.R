#' @keywords internal
# Evaluate loglikelihood for intercept only MNL, bassed on P. Rossi 2004
llmnl_int <- function(int, Obs, n_cat) {
  n_Obs 		<- length(Obs)
  Xint 		  <- matrix(c(0, int), byrow = T, ncol = n_cat, nrow = n_Obs)
  ind			  <- cbind(c(1:n_Obs), Obs)
  Xby			  <- Xint[ind]
  Xint		  <- exp(Xint)
  iota		  <- c(rep(1, (n_cat)))
  denom		  <- log(Xint %*% iota)
  return(sum(Xby - denom))
}

#' @keywords internal
# Obtain fractional log likelihood for multinomial intercept only model, bassed on P. Rossi 2004
llmnl_int_frac <- function(int, Obs, n_cat, pooled_likel, w, wgt){
  return((1 - w) * llmnl_int(int = int, Obs = Obs, n_cat = n_cat) + w * wgt * pooled_likel)
}
