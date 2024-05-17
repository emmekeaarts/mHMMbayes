#' @keywords internal
# Whenever you use C++ code in your package, you need to clean up after yourself
# when your package is unloaded. This function unloads the DLL (H. Wickham(2019). R packages)
.onUnload <- function (libpath) {
  library.dynam.unload("mHMMbayes", libpath)
}

#' @keywords internal
# simple functions used in mHMM
dif_matrix <- function(rows, cols, data = NA){
  return(matrix(data, ncol = cols, nrow = rows))
}

#' @keywords internal
nested_list <- function(n_dep, m){
  return(rep(list(vector("list", n_dep)),m))
}

#' @keywords internal
dif_vector <- function(x){
  return(numeric(x))
}

#' @keywords internal
is.whole <- function(x) {
  return(is.numeric(x) && floor(x) == x)
}

#' @keywords internal
is.mHMM <- function(x) {
  inherits(x, "mHMM")
}

#' @keywords internal
is.mHMM_tv <- function(x) {
  inherits(x, "mHMM_tv")
}

#' @keywords internal
is.mHMM_gamma <- function(x) {
  inherits(x, "mHMM_gamma")
}

#' @keywords internal
is.mHMM_prior_gamma <- function(x) {
  inherits(x, "mHMM_prior_gamma")
}

#' @keywords internal
is.mHMM_prior_emiss <- function(x) {
  inherits(x, "mHMM_prior_emiss")
}

#' @keywords internal
is.cat <- function(x) {
  inherits(x, "cat")
}

#' @keywords internal
is.cont <- function(x) {
  inherits(x, "cont")
}

#' @keywords internal
is.count <- function(x) {
  inherits(x, "count")
}

#' @keywords internal
is.mHMM_pdRW_gamma <- function(x) {
  inherits(x, "mHMM_pdRW_gamma")
}

#' @keywords internal
is.mHMM_pdRW_emiss <- function(x) {
  inherits(x, "mHMM_pdRW_emiss")
}

#' @keywords internal
hms <- function(t){
  paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"),
        formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
        formatC(t %% 60, width = 2, format = "d", flag = "0"),
        sep = ":")
}

#' @keywords internal
depth <- function(x,xdepth=0){
  if(!is.list(x)){
    return(xdepth)
  }else{
    return(max(unlist(lapply(x,depth,xdepth=xdepth+1))))
  }
}

#' @keywords internal
# Calculates the between subject variance from logmu and logvar:
logvar_to_var <- function(logmu, logvar){
  abs(exp(logvar)-1)*exp(2*logmu+logvar)
}
