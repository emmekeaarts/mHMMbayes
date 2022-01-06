#' @keywords internal
# obtain probabilities of gamma as a function of time varying covariates
# returns gamma as a matrix with elements of gamma in the columns and time in rows.
# Hence, m * m columns, and n_t rows.

# int is a list with number of elements equal to m,
# each element contains a vector with the intercepts corresponding to the i^th row in the
# transition probability matrix
# bet is a list with number of elements equal to m,
# each element contains a vector with the beta's relating to the time varying covariates,
# where the i^th element corresponds to the i^th row in the transition probability matrix
# xx_t is a vector with the time varying covariates
# m denotes the number of states

timevar_gamma <- function(int, bet, xx_t, m){

  n_xx_t <- length(xx_t)
  xx_t_int <- cbind(rep(1, n_xx_t), xx_t)
  out <- matrix(, nrow = n_xx_t, ncol = m*m)

  for(i in 1:m){
    int_bet <- c(int[[i]], bet[[i]])
    Imod <- matrix(0, m, m - 1)
    Imod[-1,] <-  diag(m - 1)
    X <- xx_t_int %x% Imod
    k <- ncol(X)
    Xbeta <- X %*% t(int_bet)
    Xbeta <- matrix(Xbeta, byrow = T, ncol = m)
    Xbeta <- exp(Xbeta)
    iota <- c(rep(1, m))
    denom <- Xbeta %*% iota
    Prob <- Xbeta / as.vector(denom)
    out[, (1 + (i-1) * m): (m + (i-1) * m)] <- Prob
  }
  return(out)
}



