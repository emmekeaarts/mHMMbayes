#' @keywords internal
# one run of the random walk metropolis sampler for an intercept only multinomial distribution
# this means no covariates at the lower/time level

mnl_RW_once <- function(int1, Obs, n_cat, mu_int_bar1, V_int1, scalar, candcov1) {
  # obtain likelihood and transition prob with the parameters sampled in the previous iteration and current sampled state sequence
  oldloglike	 		<- llmnl_int(int = int1, Obs = Obs, n_cat = n_cat)
  oldpostlike	 		<- oldloglike + dmvnorm(int1, mu_int_bar1, V_int1, log = TRUE)
  probold				  <- int_to_prob(matrix(int1, nrow = 1))

  # obtain new parameters for gamma from proposal distribution plus new likelihood
  int_new		 		  <- int1 + rmvnorm(1, rep(0, (n_cat - 1)), scalar^2 * candcov1, method = "svd")
  newloglike	 		<- llmnl_int(int = int_new, Obs = Obs, n_cat = n_cat)
  newpostlike	 		<- newloglike + dmvnorm(int_new, mu_int_bar1, V_int1, log = TRUE)
  probnew				  <- int_to_prob(matrix(int_new, nrow = 1))

  # determine to use the updated or current (previous iteration) gamma values of the parameters
  acc 				   <- min(log(1), (newpostlike - oldpostlike))
  if(acc < log(1)) {
    unif         <- log(runif(1))
  } else {
    unif         <- log(1)
  }
  if (unif <= acc) {
    draw_int		<- int_new
    accept			<- 1
    prob			  <- probnew
  } else {
    draw_int		<- int1
    accept			<- 0
    prob			  <- probold
  }
  return(list(draw_int = draw_int, accept = accept, prob = prob))
}
