#' @keywords internal
#' one run of the random walk metropolis sampler to draw an individual-, state-
#' specific poisson mean from a logNormal prior. Note that the algorithm was
#' implemented using dnorm() to sample log-means which in turn must be
#' exponentiated.

pois_RW_once <- function(lambda, Obs, mu_bar1, V_1, scalar, candcov1) {

  # obtain likelihood and transition prob with the parameters sampled in the previous iteration and current sampled state sequence
  oldloglike	 	<- llpois(lambda = lambda, Obs = Obs)
  oldpostlike	 	<- oldloglike + dnorm(log(lambda), mu_bar1, V_1, log = TRUE)

  # obtain new parameters for gamma from proposal distribution plus new likelihood
  lambda_new		<- exp(log(lambda) + rnorm(1, 0, scalar * sqrt(candcov1)))
  newloglike	 	<- llpois(lambda = lambda_new, Obs = Obs)
  newpostlike	 	<- newloglike + dnorm(log(lambda_new), mu_bar1, V_1, log = TRUE)

  # determine to use the updated or current (previous iteration) gamma values of the parameters
  acc 			<- min(log(1), (newpostlike - oldpostlike))
  if(acc < log(1)) {
    unif        <- log(runif(1))
  } else {
    unif        <- log(1)
  }
  if (unif <= acc) {
    draw_lambda	<- lambda_new
    accept		<- 1
  } else {
    draw_lambda	<- lambda
    accept		<- 0
  }

  return(list(draw_lambda = draw_lambda, accept = accept))
}
