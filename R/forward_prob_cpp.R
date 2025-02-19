#' @keywords internal
# Calculates the probabilities of observing each state at each point in time given
# the observations of all dependent variables, used for the forward probabilities
# Based on Zuchini 2016.
all1 <- function(x, emiss, n_dep, data_distr){
  inp <- rep(list(NULL), n_dep)
  if(data_distr == "categorical"){
    for(q in 1:n_dep){
      inp[[q]] <- t(emiss[[q]][,x[,q]])
    }
  } else if (data_distr == "continuous"){
    for(q in 1:n_dep){
      inp[[q]] <- outer(x[,q], Y = emiss[[q]][,1], FUN = stats::dnorm, sd = rep(sqrt(emiss[[q]][,2]), each = dim(x)[1]))
    }
  } else if (data_distr == "count") {
    for(q in 1:n_dep){
      inp[[q]] <- outer(x[,q], emiss[[q]][,1], FUN = stats::dpois)
    }
  }
  allprobs <- Reduce("*", inp)
  return(allprobs)
}


#' @keywords internal
# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.
cat_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, gamma, delta = NULL){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "categorical")
  out <- cat_mult_fw_cpp(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
  return(out)
}

#' @keywords internal
# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.
cont_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, gamma, delta = NULL){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "continuous")
  out <- cat_mult_fw_cpp(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
  return(out)
}

#' @keywords internal
# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.
count_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, gamma, delta = NULL){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "count")
  out <- cat_mult_fw_cpp(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
  return(out)
}



#' Calculates the forward probabilities, used for sampling the state sequence
#' specific for when data includes a time varying covariate
#'
#' @export
#'

cat_tv_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, tgamma, delta = NULL, data_distr){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - matrix(tgamma[1,], byrow = TRUE, ncol = m) + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "categorical")
  out <- cat_tv_mult_fw_cpp(allprobs = allprobs, tgamma = tgamma, m = m, n = n, delta = delta)
  return(out)
}


#' Calculates the forward probabilities, used for sampling the state sequence
#' specific for when data includes a time varying covariate
#'
#' @export
#'

cont_tv_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, tgamma, delta = NULL, data_distr){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - matrix(tgamma[1,], byrow = TRUE, ncol = m) + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "continuous")
  out <- cat_tv_mult_fw_cpp(allprobs = allprobs, tgamma = tgamma, m = m, n = n, delta = delta)
  return(out)
}


#' Calculates the forward probabilities, used for sampling the state sequence
#' specific for when data includes a time varying covariate
#'
#' @export
#'

count_tv_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, tgamma, delta = NULL, data_distr){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - matrix(tgamma[1,], byrow = TRUE, ncol = m) + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "count")
  out <- cat_tv_mult_fw_cpp(allprobs = allprobs, tgamma = tgamma, m = m, n = n, delta = delta)
  return(out)
}
