#' @keywords internal
# computes probabilities from intercepts
int_to_prob <- function(int1) {
  exp_int1 	<- matrix(exp(c(0, int1)), nrow  = 1)
  prob1 		<- exp_int1 / as.vector(exp_int1 %*% c(rep(1, (dim(exp_int1)[2]))))
  return(prob1)
}
