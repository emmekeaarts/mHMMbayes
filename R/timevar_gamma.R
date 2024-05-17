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
# xx_t is a matrix, with the first column a vector of 1's denoting the intercept, and
# the remaining column the time varying covariate(s)
# m denotes the number of states

timevar_gamma <- function(int, bet, xx_t, m){

  n_xx_t <- dim(xx_t)[1]
  out <- matrix(NA_real_, nrow = n_xx_t, ncol = m*m)

  for(i in 1:m){
    # int_bet <- matrix(c(int[i,], bet[i,]), ncol = 1)
    int_bet <- matrix(c(int[i,], as.numeric(bet[[i]])), ncol = 1)
    Imod <- matrix(0, m, m - 1)
    Imod[-1,] <-  diag(m - 1)
    X <- xx_t %x% Imod
    k <- ncol(X)
    Xbeta <- X %*% int_bet
    Xbeta <- matrix(Xbeta, byrow = T, ncol = m)
    Xbeta <- exp(Xbeta)
    iota <- c(rep(1, m))
    denom <- Xbeta %*% iota
    Prob <- Xbeta / as.vector(denom)
    out[, (1 + (i-1) * m): (m + (i-1) * m)] <- Prob
  }
  return(out)
}
