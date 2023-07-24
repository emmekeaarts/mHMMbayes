#' Summary of a \code{mHMM} object
#'
#' @param object A \code{mHMM} object.
#' @param ... Other parameters passed down to \code{summary()}.
#'
#' @export

summary.mHMM <- function(object, ...){
  input   <- object$input
  dep_labels <- input$dep_labels
  n_subj  <- input$n_subj
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  q_emiss <- input$q_emiss
  n_dep   <- input$n_dep
  type_covar <- input$covar_type
  gamma_int <- matrix(apply(object$gamma_int_bar[((burn_in + 1): J),], 2, median),
                      byrow = TRUE, ncol = m-1, nrow = m)
  gamma_pop <- round(int_to_prob(gamma_int),3)
  colnames(gamma_pop) <- paste("To state", 1:m)
  rownames(gamma_pop) <- paste("From state", 1:m)
  cat("State transition probability matrix","\n",  "(at the group level):", "\n", "\n")
  print(gamma_pop)
  cat("\n")
  if (!is.null(type_covar[[1]])){
    regr_coeff_gamma <- data.frame(From_state = c(rep(1:m, each = m - 1)),
                                   To_state = rep(2:m, m))|>
      dplyr::mutate(Beta = round(apply(object$gamma_cov_bar[(burn_in + 1): J,],
                                       2, median),3),
                    CrI_lower = round(apply(object$gamma_cov_bar[(burn_in + 1): J,],
                                            2, quantile, probs = 0.025),3),
                    CrI_upper = round(apply(object$gamma_cov_bar[(burn_in + 1): J,],
                                            2, quantile, probs = 0.975),3),
                    ` ` = ifelse(CrI_lower < 0 & CrI_upper > 0, " ", "*"))
    cat("Regression coefficients predicting the transition probabilities","\n",
        "(at the group level):", "\n", "\n")
    print(regr_coeff_gamma, row.names = F)
    cat("Note: [*] 95% credible interval does not include zero.")
    cat("\n", "\n")
  }
  cat("Emission distribution for each of the dependent variables","\n",
      "(at the group level):", "\n", "\n")
  EM_int <- EM_pop <- vector("list", n_dep)
  names(EM_pop) <- dep_labels
  for(i in 1:n_dep){
    EM_int[[i]] <- matrix(apply(object$emiss_int_bar[[i]][((burn_in + 1): J),],
                                2, median), byrow = TRUE, ncol = q_emiss[i]-1, nrow = m)
    EM_pop[[i]] <- round(int_to_prob(EM_int[[i]]),3)
    colnames(EM_pop[[i]]) <- paste("Category", 1:q_emiss[i])
    rownames(EM_pop[[i]]) <- paste("State", 1:m)
  }
  print(EM_pop)
  cat("\n")
  # assume covar lists for emiss are the same
  if (!is.null(type_covar[[2]])){
    regr_coeff_emiss <- vector("list", n_dep)
    names(regr_coeff_emiss) <- dep_labels
    for (i in 1:n_dep){
      regr_coeff_emiss[[i]] <- data.frame(Category = c(rep(2:q_emiss[i], m)),
                                          State = rep(1:m, each = q_emiss[i]-1))|>
        dplyr::mutate(Beta = round(apply(object$emiss_cov_bar[[i]][(burn_in + 1): J,],
                                         2, median),3),
          CrI_lower = round(apply(object$emiss_cov_bar[[i]][(burn_in + 1): J,],
                                  2, quantile, probs = 0.025),3),
          CrI_upper = round(apply(object$emiss_cov_bar[[i]][(burn_in + 1): J,],
                                  2, quantile, probs = 0.975),3),
          ` ` = ifelse(CrI_lower < 0 & CrI_upper > 0, " ", "*")
          )
    }
    cat("Regression coefficients predicting the emission probabilities for each of the dependent variables","\n",  "(at the group level):", "\n", "\n")
    print(regr_coeff_emiss, row.names = F)
    cat("Note: [*] 95% credible interval does not include zero.")
  }
  cat("\n")
}


utils::globalVariables(c("quantile", "CrI_lower", "CrI_upper", "Beta", "To_state", "State"))
