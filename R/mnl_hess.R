# Obtain mnl -Expected[Hessian]  for intercept only model, bassed on P.Rossi 2004
#' @keywords internal
mnlHess_int <- function(int, Obs, n_cat){
  n_Obs 	<- length(Obs)
  betas   <- matrix(c(0, int), byrow = T, ncol = n_cat)
  prob    <- exp(betas) / sum(exp(betas))
  Hess    <- (diag(x = prob[-1], nrow = n_cat-1) - outer(prob[-1],prob[-1])) * n_Obs
  return(Hess)
}

#' @keywords internal
# Obtain mnl -Expected[Hessian]  for model including time variant cov, bassed on P.Rossi 2004
mnlHess_bet <- function(bet, Obs, xx_t, n_cat){
  n_Obs <- length(Obs)
  Imod <- matrix(0, n_cat, n_cat - 1)
  Imod[-1,] <-  diag(n_cat - 1)
  X <- xx_t %x% Imod
  k <- ncol(X)
  Xbeta <- X %*% bet
  Xbeta <- matrix(Xbeta, byrow = T, ncol = n_cat)
  Xbeta <- exp(Xbeta)
  iota <- c(rep(1, n_cat))
  denom <- Xbeta %*% iota
  Prob <- Xbeta / as.vector(denom)
  Hess    <- matrix(double(k * k), ncol = k)

  for (hi in 1:n_Obs) {
    p     <- as.vector(Prob[hi, ])
    A     <- diag(p) - outer(p, p)
    Xt    <- X[(n_cat*(hi - 1) + 1):(n_cat * hi),]
    Hess=Hess+crossprod(Xt,A)%*%Xt

    Xt    <- diag(n_cat)[, -1]
    Hess  <- Hess + crossprod(Xt, A) %*% Xt
  }
  return(Hess)
}
