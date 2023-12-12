#' @keywords internal
# Calculates the between subject variance from logmu and logvar:
obtain_varmu <- function(logmu, logvar){
  abs(exp(logvar)-1)*exp(2*logmu+logvar)
}

#' @keywords internal
# Calculates necessary logvar (log scale) to get the desired between subject
# variance in the real domain (depends also on the emission means chosen):
obtain_logvar <- function(mu, varmu){
  logmu = log(mu)
  log(0.5*exp(-2*logmu)*(exp(2*logmu) + sqrt(4*exp(2*logmu)*varmu+exp(4*logmu))))
}
