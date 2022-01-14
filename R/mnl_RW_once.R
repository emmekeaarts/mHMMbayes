#' @keywords internal
# one run of the random walk metropolis sampler for an intercept only multinomial distribution
# this means no covariates at the lower/time level

mnl_RW_once <- function(int1, Obs, n_cat, mu_int_bar1, V_int1, scalar, candcov1) {
  # obtain likelihood and transition prob with the parameters sampled in the previous iteration and current sampled state sequence
  oldloglike	 		<- llmnl_int(int = int1, Obs = Obs, n_cat = n_cat)
  oldpostlike	 		<- oldloglike + dmvnorm(int1, mu_int_bar1, V_int1, log = TRUE)
  probold				  <- int_to_prob(int1)

  # obtain new parameters for gamma from proposal distribution plus new likelihood
  int_new		 		  <- int1 + rmvnorm(1, rep(0, (n_cat - 1)), scalar^2 * candcov1, method = "svd")
  newloglike	 		<- llmnl_int(int = int_new, Obs = Obs, n_cat = n_cat)
  newpostlike	 		<- newloglike + dmvnorm(int_new, mu_int_bar1, V_int1, log = TRUE)
  probnew				  <- int_to_prob(int_new)

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


#' @keywords internal
# one run of the random walk metropolis sampler for an multinomial distribution
# with a random intercept (estimated here), and a fixed beta for a time-varying covariate,
# which is not estimated at this step but needs to be incorporated when obtaining the
# log likelihood estimates

mnl_RW_once_int_fixed_bet <- function(int1, Obs, xx_t, n_cat, mu_int_bar1, bet, V_int1, scalar, candcov1) {
  # obtain likelihood and transition prob with the parameters sampled in the previous iteration and current sampled state sequence
  oldloglike	 		<- llmnl_bet(bet = c(int1, bet), Obs = Obs, xx_t = xx_t, n_cat = n_cat)
  oldpostlike	 		<- oldloglike + dmvnorm(int1, mu_int_bar1, V_int1, log = TRUE)
  probold				  <- int_to_prob(int1)

  # obtain new parameters for gamma from proposal distribution plus new likelihood
  int_new		 		  <- int1 + rmvnorm(1, rep(0, (n_cat - 1)), scalar^2 * candcov1, method = "svd")
  newloglike	 		<- llmnl_bet(bet = c(int_new, bet), Obs = Obs, xx_t = xx_t, n_cat = n_cat)
  newpostlike	 		<- newloglike + dmvnorm(int_new, mu_int_bar1, V_int1, log = TRUE)
  probnew				  <- int_to_prob(int_new)

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


#' @keywords internal
# one run of the random walk metropolis sampler for an multinomial distribution
# with a fixed beta for (one) time-varying covariate, in addition to random intercept.
# The random intercepts are not estimated here, but incorporated when obtaining the
# log likelihood estimates for the fixed beta

mnl_RW_once_fixed_bet_pooled <- function(bet1, Obs_p, xx_t_p, n_cat, mu_bet_bar1, int_p, V_bet1, scalar, candcov2, n_subj) {
  # obtain likelihood and transition prob with the parameters sampled in the previous iteration and current sampled state sequence
  bet_p           <- cbind(int_p, matrix(bet1, byrow = TRUE, ncol = n_cat - 1, nrow = n_subj))
  oldloglike	 		<- llmnl_bet_pooled(bet_p = bet_p, Obs_p = Obs_p, xx_t_p = xx_t_p, n_cat = n_cat, n_subj = n_subj)
  oldpostlike	 		<- oldloglike + dmvnorm(unlist(bet1), mu_bet_bar1, V_bet1, log = TRUE)

  # obtain new parameters for gamma from proposal distribution plus new likelihood
  bet_new		 		  <- bet1 + rmvnorm(1, rep(0, (n_cat - 1)), scalar^2 * candcov2, method = "svd")
  bet_p_new       <- cbind(int_p, matrix(bet_new, byrow = TRUE, ncol = n_cat - 1, nrow = n_subj))
  newloglike	 		<- llmnl_bet_pooled(bet_p = bet_p_new, Obs_p = Obs_p, xx_t_p = xx_t_p, n_cat = n_cat, n_subj = n_subj)
  newpostlike	 		<- newloglike + dmvnorm(bet_new, mu_bet_bar1, V_bet1, log = TRUE)

  # determine to use the updated or current (previous iteration) gamma values of the parameters
  acc 				   <- min(log(1), (newpostlike - oldpostlike))
  if(acc < log(1)) {
    unif         <- log(runif(1))
  } else {
    unif         <- log(1)
  }
  if (unif <= acc) {
    draw_int		<- bet_new
    accept			<- 1
  } else {
    draw_int		<- bet1
    accept			<- 0
  }
  return(list(draw_int = draw_int, accept = accept))
}
