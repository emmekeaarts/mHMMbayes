#' @export
#'
print.mHMM <- function(x, ...){
  if(sum(objects(x$PD_subj[[1]]) %in% "log_likl") != 1){
    stop("The input object is created using an earlier version of the mHMMbayes package. Please re-run the function mHMM with the current package version, or post-process the object using the earlier version of the package.")
  }
  input   <- x$input
  n_subj  <- input$n_subj
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  n_vary  <- input$n_vary
  data_distr <- input$data_distr
  n_dep   <- input$n_dep
  cat("Number of subjects:", n_subj, "\n")
  cat("\n")
  cat(J, "iterations used in the MCMC algorithm with a burn in of", burn_in, "\n")
  LL      <- numeric(n_subj)
  for(s in 1:n_subj){
    LL[s] <- median(x$PD_subj[[s]]$log_likl[((burn_in + 1): J), 1])
  }
  if(data_distr == 'categorical'){
    q_emiss <- input$q_emiss
    n_par <- sum((q_emiss-1)*m)+(m-1)*m
    AIC   <- 2 * n_par - (2 * LL)
    AICc  <- ((2 * n_vary * n_par) / (n_vary - n_par - 1)) - (2 * LL)
  } else if (data_distr == 'continuous'){
    n_par <- m * n_dep * 2 + (m - 1) * m
    AIC   <- 2 * n_par - (2 * LL)
    AICc   <- ((2 * n_vary * n_par) / (n_vary - n_par - 1)) - (2 * LL)
  } else if (data_distr == 'count'){
    n_par <- m * n_dep + (m - 1) * m
    AIC   <- 2 * n_par - (2 * LL)
    AICc  <- ((2 * n_vary * n_par) / (n_vary - n_par - 1)) - (2 * LL)
  }

  cat("Average Log likelihood over all subjects:", mean(LL), "\n")
  cat("Average AIC over all subjects: ", mean(AIC), "\n")
  cat("Average AICc over all subjects:", mean(AICc), "\n")
  cat("\n")
  cat("Number of states used:", m, "\n")
  cat("\n")
  cat("Number of dependent variables used:", n_dep, "\n")
  cat("\n")
  cat("Type of dependent variable(s):", data_distr, "\n")
  cat("\n")
}
