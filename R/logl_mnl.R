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

#' @keywords internal
# Obtain xxx, bassed on P. Rossi 2004
## xx_t is matrix with observed time varying covariates, with first colom being a vector of 1's denoting the intercept
## note that for xx_t, take care that only the observed output is given for time points where
## we are using the ll calculation for, so for gamma, only for points at the sampled state sequence where state == i
## so xx_t has to be the same lenght as Obs

llmnl_bet <- function(bet, Obs, xx_t, n_cat){
  n_Obs <- length(Obs)
  Imod <- matrix(0, n_cat, n_cat - 1)
  Imod[-1,] <-  diag(n_cat - 1)
  X <- xx_t %x% Imod
  Xbeta <- X %*% bet
  Xbeta <- matrix(Xbeta, byrow = T, ncol = n_cat)
  ind <- cbind(c(1:n_Obs), Obs)
  xby <- Xbeta[ind]
  Xbeta <- exp(Xbeta)
  iota <- c(rep(1, n_cat))
  denom <- log(Xbeta %*% iota)
  return(sum(xby - denom))
}


#' @keywords internal
# Obtain fractional log likelihood for multinomial model with covariates, bassed on P. Rossi 2004
llmnl_bet_frac <- function(bet, Obs, xx_t, n_cat, pooled_likel, w, wgt){
  return((1 - w) * llmnl_bet(bet = bet, Obs = Obs, xx_t = xx_t, n_cat = n_cat) + w * wgt * pooled_likel)
}


#' @keywords internal
# Obtain xxx, bassed on P. Rossi 2004
## xx_t is matrix with observed time varying covariates, with first colom being a vector of 1's denoting the intercept
## note that for xx_t, take care that only the observed output is given for time points where
## we are using the ll calculation for, so for gamma, only for points at the sampled state sequence where state == i
## so xx_t has to be the same lenght as Obs2

## bet2 is a matrix, with in the first n_cat - 1 columns the intercepts that vary over the subjects,
## and in the remaining n_cat - 1 colums (currently we only support incorporating ONE time varying covariant)
## is the current estimate of the beta's which are fixed over the subjects
## Obs2 is a matrix, with subject indicator in the first colomn, and the observed outcomes in the second column
llmnl_bet_pooled <- function(bet_p, Obs_p, xx_t_p, n_cat, n_subj){
  n_Obs <- nrow(Obs_p)
  Imod <- matrix(0, n_cat, n_cat - 1)
  Imod[-1,] <-  diag(n_cat - 1)
  Xbeta <- NULL
  for(i in 1:n_subj){
    X <- xx_t_p[Obs_p[,1] == i,] %x% Imod
    Xbeta <- rbind(Xbeta, X %*% bet_p[i,])
  }
  Xbeta <- matrix(Xbeta, byrow = T, ncol = n_cat)
  ind <- cbind(c(1:n_Obs), Obs_p[,2])
  xby <- Xbeta[ind]
  Xbeta <- exp(Xbeta)
  iota <- c(rep(1, n_cat))
  denom <- log(Xbeta %*% iota)
  return(sum(xby - denom))
}
