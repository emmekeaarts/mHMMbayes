#' @keywords internal
# computes probabilities from intercepts
int_to_prob <- function(int1) {
  if(is.matrix(int1)){
    prob1 <- matrix(nrow = nrow(int1), ncol = ncol(int1) + 1)
    for(r in 1:nrow(int1)){
      exp_int1 	<- matrix(exp(c(0, int1[r,])), nrow  = 1)
      prob1[r,] <- exp_int1 / as.vector(exp_int1 %*% c(rep(1, (dim(exp_int1)[2]))))
    }
  } else {
    exp_int1 	<- matrix(exp(c(0, int1)), nrow  = 1)
    prob1 		<- exp_int1 / as.vector(exp_int1 %*% c(rep(1, (dim(exp_int1)[2]))))
  }
  return(round(prob1,4))
}

# computes intercepts from probabilities, per row of input matrix
# first catagory is reference catagory
prob_to_int <- function(prob1){
  prob1 <- prob1 + 0.00001
  b0 <- matrix(NA, nrow(prob1), ncol(prob1)-1)
  sum_exp <- numeric(nrow(prob1))
  for(r in 1:nrow(prob1)){
    sum_exp[r] <- (1/prob1[r,1]) - 1
    for(cr in 2:ncol(prob1)){
      #for every b0 except the first collumn (e.g. b012 <- log(y12/y11-y12))
      b0[r,(cr-1)] <- log(prob1[r,cr]*(1+sum_exp[r]))
    }
  }
  return(round(b0,4))
}
