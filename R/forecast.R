#' Forecasts (forward predictions) using a multilevel hidden Markov model
#'
#' Ideas:
#' (1) use MAPs similarly as done with the Viterbi algorithm:
#'  +++ fast
#'  --- no measure of uncertainty on forcasts
#' (2) obtain forecasts over iterations (after burn-in):
#'  +++ uncertainty on forecasts
#'  --- slow
#'
#' To add:
#'  - option to request simulating (forecasting) observations or at least
#'    predicting the probability of observations for a number of timesteps.
#'
#' @examples
#'
#'
#' @export

forecast_mHMM1 <- function(s_data, object, forecast_steps = 1, burn_in = NULL) {

  if (!is.mHMM(object)){
    stop("The input object used should be from the class mHMM, obtained by using the function mHMM.")
  }
  if(sum(objects(object$PD_subj[[1]]) %in% "log_likl") != 1){
    stop("The input object is created using an earlier version of the mHMMbayes package. Please re-run the function mHMM with the current package version, or post-process the object using the earlier version of the package.")
  }
  input      <- object$input
  data_distr <- input$data_distr
  id         <- unique(s_data[,1])
  n_subj     <- length(id)
  if(length(object$PD_subj) != n_subj){
    stop("s_data used should be from the same subjects used for creating the object in mHMM.
         The number of subjects in the datasets are not the same.")
  }
  n_vary     <- table(s_data[,1])
  max_n      <- max(n_vary)
  state_seq  <- matrix(,ncol = n_subj, nrow = max_n)
  probs      <- vector(mode = "list", length = n_subj)
  n_dep      <- input$n_dep
  m          <- input$m
  if(is.null(burn_in)){
    burn_in  <- input$burn_in
  }
  J          <- input$J
  if (burn_in >= (J-1)){
    stop(paste("The specified burn in period should be at least 2 points smaller
               compared to the number of iterations J, J =", J))
  }

  # Get subject-specific emissions
  if(data_distr == "categorical"){
    q_emiss    <- input$q_emiss
    int_est_emiss <- rep(list(lapply(q_emiss-1, dif_matrix, rows = m)), n_subj)
    est_emiss  <- rep(list(lapply(q_emiss, dif_matrix, rows = m)), n_subj)
    for(s in 1:n_subj){
      for(j in 1:n_dep){
        int_est_emiss[[s]][[j]][] <- matrix(apply(object$emiss_int_subj[[s]][[j]][burn_in:J, ], 2, median),
                                            byrow = TRUE, ncol = q_emiss[j]-1, nrow = m)
        est_emiss[[s]][[j]][] <- int_to_prob(int_est_emiss[[s]][[j]])
      }
    }
  } else if(data_distr == "continuous"){
    est_emiss  <- rep(list(rep(list(matrix(NA_real_,nrow = m, ncol = 2)),n_dep)), n_subj)
    for(s in 1:n_subj){
      for(q in 1:n_dep){
        est_emiss[[s]][[q]][] <- matrix(round(c(apply(object$PD_subj[[s]]$cont_emiss[((burn_in + 1): J),((q-1) * m + 1):(q * m)], 2, median),
                                                apply(object$PD_subj[[s]]$cont_emiss[((burn_in + 1): J), (n_dep * m + (q-1) * m + 1):(n_dep * m + q * m)], 2, median)),3),
                                        ncol = 2, nrow = m)
      }
    }
  } else if(data_distr == "count"){
    est_emiss  <- rep(list(rep(list(matrix(NA_real_,nrow = m, ncol = 1)),n_dep)), n_subj)
    for(s in 1:n_subj){
      for(q in 1:n_dep){
        est_emiss[[s]][[q]][] <- matrix(round(apply(object$PD_subj[[s]]$count_emiss[((burn_in + 1): J),((q-1) * m + 1):(q * m)], 2, median),3),
                                        ncol = 1, nrow = m)
      }
    }
  }

  # Get subject-specific transitions
  est_gamma <- obtain_gamma(object, level = "subject")

  # Obtain the forward probabilities and make forecasts

  for(s in 1:n_subj){
    emiss   <- est_emiss[[s]]
    gamma   <- est_gamma[[s]]
    if(data_distr == "categorical"){
      probs[[s]]    <- cat_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                                      matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
                                            m = m, emiss = emiss, gamma = gamma, n_dep = n_dep, delta=NULL)[[1]][,(n_vary[s]+1):(n_vary[s]+forecast_steps)]
    } else if(data_distr == "continuous"){
      probs[[s]]    <- cont_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                                       matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
                                             m = m, emiss = emiss, n_dep = n_dep, gamma = gamma)[[1]][,(n_vary[s]+1):(n_vary[s]+forecast_steps)]
    } else if(data_distr == "count"){
      probs[[s]]    <- count_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                                        matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
                                              m = m, emiss = emiss, n_dep = n_dep, gamma = gamma)[[1]][,(n_vary[s]+1):(n_vary[s]+forecast_steps)]
    }

    # # Sample states
    # state_seq[1:n_vary[s], s] <- apply(probs[[s]], 2, which.max)

    # # Using the forward probabilites, sample the state sequence in a backward manner.
    # # In addition, saves state transitions in trans, and conditional observations within states in cond_y
    # # trans[[s]]					                  <- vector("list", m)
    # sample_path[[s]][n_vary[[s]], iter] 	<- sample(1:m, 1, prob = c(alpha[, n_vary[[s]]]))
    # for(t in (subj_data[[s]]$n_t - 1):1){
    #   sample_path[[s]][t,iter] 	              <- sample(1:m, 1, prob = (alpha[, t] * gamma[[s]][,sample_path[[s]][t + 1, iter]]))
    #   # trans[[s]][[sample_path[[s]][t,iter]]]	<- c(trans[[s]][[sample_path[[s]][t, iter]]], sample_path[[s]][t + 1, iter])
    # }

    # for(h in 1:forecast_steps){
    #
    # }

  }

  # Return results
  return(forecast_probs = probs)

}

# # Example usage:
# set.seed(123)
# # Simulate data using the mHMMbayes package
# num_states <- 2
# num_obs <- 5
# sim_data <- sim_mHMM(NStates = num_states, n = num_obs, nchains = 3)
#
# # Fit a multilevel Bayesian HMM using mHMMbayes
# fit <- fit_mHMM(Y = sim_data$Y, nchains = 3)
#
# # Perform forward prediction
# num_steps <- 3
# result <- forward_prediction_multilevel(fit, num_steps)
#
# # Print the results
# cat("True States:", sim_data$state[1, , ], "\n")
# cat("Observed Data:", sim_data$Y[1, , ], "\n")
# cat("Predicted States:\n", result$predicted_states[1, , , drop = FALSE], "\n")
# cat("Predicted Observations:\n", result$predicted_observations[1, , drop = FALSE], "\n")

### Example on continuous simulated data
library(mHMMbayes)

n_t     <- 500
n       <- 10
m       <- 3
n_dep   <- 2

gamma   <- matrix(c(0.99, 0.005, 0.005,
                    0.08, 0.9, 0.02,
                    0.05, 0.15, 0.8), ncol = m, byrow = TRUE)

emiss_distr <- list(matrix(c( 50, 10,
                              100, 10,
                              150, 10), nrow = m, byrow = TRUE),
                    matrix(c(5, 2,
                             10, 5,
                             20, 3), nrow = m, byrow = TRUE))

data_cont <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous', gen = list(m = m, n_dep = n_dep),
                  gamma = gamma, emiss_distr = emiss_distr, var_gamma = .1, var_emiss = c(.5, 0.01))

# Specify hyper-prior for the continuous emission distribution
manual_prior_emiss <- prior_emiss_cont(
                        gen = list(m = m, n_dep = n_dep),
                        emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
                                         matrix(c(7, 8, 18), nrow = 1)),
                        emiss_K0 = list(1, 1),
                        emiss_V =  list(rep(100, m), rep(25, m)),
                        emiss_nu = list(1, 1),
                        emiss_a0 = list(rep(1, m), rep(1, m)),
                        emiss_b0 = list(rep(1, m), rep(1, m)))

# Run the model on the simulated data:
# Note that for reasons of running time, J is set at a ridiculous low value.
# One would typically use a number of iterations J of at least 1000,
# and a burn_in of 200.
out_3st_cont_sim <- mHMM(s_data = data_cont$obs,
                         data_distr = 'continuous',
                         gen = list(m = m, n_dep = n_dep),
                         start_val = c(list(gamma), emiss_distr),
                         emiss_hyp_prior = manual_prior_emiss,
                         mcmc = list(J = 500, burn_in = 250))

summary(out_3st_cont_sim)

t(mHMMbayes::forecast_mHMM1(object = out_3st_cont_sim, s_data = data_cont$obs, forecast_steps = 50)[[3]])
