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
is.mHMM_cont <- function(x) {
  inherits(x, "mHMM_cont")
}

#' @keywords internal
is.mHMM_vary <- function(x) {
  inherits(x, "mHMM_vary")
}

#' @keywords internal
is.mHMM_gamma <- function(x) {
  inherits(x, "mHMM_gamma")
}

#' @keywords internal
is.mHMM_emiss <- function(x) {
  inherits(x, "mHMM_emiss")
}

#' @keywords internal
is.mHMM_emiss_cont <- function(x) {
  inherits(x, "mHMM_emiss_cont")
}

#' @keywords internal
hms <- function(t){
  paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"),
        formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
        formatC(t %% 60, width = 2, format = "d", flag = "0"),
        sep = ":")
}
#' @keywords internal
isNested <- function(x) {
  out <- FALSE
  strout <- utils::capture.output(str(x))
  idx <- grep("\\$.*List", strout)
  if (length(idx)) {
    out <- TRUE
  }
  return(out)
}

#' @keywords internal
t_col <- function(color, percent = 60, name = NULL) {

  rgb.val <- col2rgb(color)

  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
    invisible(t.col)
}
