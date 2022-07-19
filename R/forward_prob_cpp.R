#' @keywords internal
# Calculates the probabilities of observing each state at each point in time given
# the observations of all dependent variables, used for the forward probabilities
# Based on Zuchini 2016.
all1 <- function(x, emiss, n_dep){
  inp <- rep(list(NULL), n_dep)
  for(q in 1:n_dep){
    inp[[q]] <- t(emiss[[q]][,x[,q]])
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
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep)
  out <- cat_mult_fw_cpp(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
  return(out)
}

