#' @keywords internal
# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.

cat_Mult_HMM_fw <- function(x, m, emiss, n_dep, gamma, delta = NULL){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  lalpha   <- alpha_prob <- matrix(NA, m, n)
  inp <- rep(list(NULL), n_dep)
  for(q in 1:n_dep){
    inp[[q]] <- t(emiss[[q]][,x[,q]])
  }
  allprobs <- Reduce("*", inp)
  foo             <- delta * allprobs[1, ]
  sumfoo          <- sum(foo)
  alpha_prob[, 1] <- foo/sumfoo
  lscale          <- log(sumfoo)
  lalpha[, 1]     <- log(alpha_prob[, 1]) + lscale
  for (i in 2:n){
    foo              <- alpha_prob[, (i - 1)] %*% gamma * allprobs[i, ]
    sumfoo           <- sum(foo)
    alpha_prob[, i]  <- foo / sumfoo
    lscale           <- lscale + log(sumfoo)
    lalpha[, i]       <- log(alpha_prob[, i]) + lscale
  }
  list(la = lalpha, forward_p = alpha_prob)
}

