#' Forecasts (forward predictions) using a multilevel hidden Markov model
#'
#' Ideas:
#' (1) use MAPs similarly as done with the Viterbi algorithm:
#'  +++ fast
#'  --- no measure of uncertainty on forecasts
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
#' ### Example on continuous simulated data
#' library(mHMMbayes)
#'
#' n_t     <- 500
#' n       <- 10
#' m       <- 3
#' n_dep   <- 2
#'
#' gamma   <- matrix(c(0.99, 0.005, 0.005,
#'                     0.08, 0.9, 0.02,
#'                     0.05, 0.15, 0.8), ncol = m, byrow = TRUE)
#'
#' emiss_distr <- list(matrix(c( 50, 10,
#'                               100, 10,
#'                               150, 10), nrow = m, byrow = TRUE),
#'                     matrix(c(5, 2,
#'                              10, 5,
#'                              20, 3), nrow = m, byrow = TRUE))
#'
#' data_cont <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous', gen = list(m = m, n_dep = n_dep),
#'                       gamma = gamma, emiss_distr = emiss_distr, var_gamma = .1, var_emiss = c(.5, 0.01))
#'
#' # Specify hyper-prior for the continuous emission distribution
#' manual_prior_emiss <- prior_emiss_cont(
#'   gen = list(m = m, n_dep = n_dep),
#'   emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
#'                    matrix(c(7, 8, 18), nrow = 1)),
#'   emiss_K0 = list(1, 1),
#'   emiss_V =  list(rep(100, m), rep(25, m)),
#'   emiss_nu = list(1, 1),
#'   emiss_a0 = list(rep(1, m), rep(1, m)),
#'   emiss_b0 = list(rep(1, m), rep(1, m)))
#'
#' # Run the model on the simulated data:
#' # Note that for reasons of running time, J is set at a ridiculous low value.
#' # One would typically use a number of iterations J of at least 1000,
#' # and a burn_in of 200.
#' out_3st_cont_sim <- mHMM(s_data = data_cont$obs,
#'                          data_distr = 'continuous',
#'                          gen = list(m = m, n_dep = n_dep),
#'                          start_val = c(list(gamma), emiss_distr),
#'                          emiss_hyp_prior = manual_prior_emiss,
#'                          mcmc = list(J = 500, burn_in = 250))
#'
#' summary(out_3st_cont_sim)
#'
#' t(forecast_mHMM1(object = out_3st_cont_sim, s_data = data_cont$obs, forecast_steps = 50)[[3]])
#'
#' @export

forecast_mHMM1 <- function(s_data, object, forecast_steps = 1, value_range = NULL, burn_in = NULL, return_all = FALSE) {

  if (!is.mHMM(object)){
    stop("The input object used should be from the class mHMM, obtained by using the function mHMM.")
  }
  if(sum(objects(object$PD_subj[[1]]) %in% "log_likl") != 1){
    stop("The input object is created using an earlier version of the mHMMbayes package.
         Please re-run the function mHMM with the current package version, or post-process the object
         using the earlier version of the package.")
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
  state_seq  <- matrix(NA, ncol = n_subj, nrow = max_n)
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

      if(return_all){
        probs[[s]]    <- t(cat_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                                          matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
                                                m = m, emiss = emiss, gamma = gamma, n_dep = n_dep, delta=NULL)[[1]])
      } else {
        probs[[s]]    <- t(cat_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                                          matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
                                                m = m, emiss = emiss, gamma = gamma, n_dep = n_dep, delta=NULL)[[1]][,(n_vary[s]+1):(n_vary[s]+forecast_steps)])
      }

      # colnames(probs[[s]]) <- paste0("pr_state_",1:m)


    } else if(data_distr == "continuous"){

      if(return_all){
        probs[[s]]    <- t(cont_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                                           matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
                                                 m = m, emiss = emiss, n_dep = n_dep, gamma = gamma)[[1]])
      } else {
        probs[[s]]    <- t(cont_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                                           matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
                                                 m = m, emiss = emiss, n_dep = n_dep, gamma = gamma)[[1]][,(n_vary[s]+1):(n_vary[s]+forecast_steps)])
      }

      # colnames(probs[[s]]) <- paste0("pr_state_",1:m)
    } else if(data_distr == "count"){

      if(return_all){
        probs[[s]]    <- t(count_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                                            matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
                                                  m = m, emiss = emiss, n_dep = n_dep, gamma = gamma)[[1]])
      } else {
        probs[[s]]    <- t(count_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                                            matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
                                                  m = m, emiss = emiss, n_dep = n_dep, gamma = gamma)[[1]][,(n_vary[s]+1):(n_vary[s]+forecast_steps)])
      }

      # colnames(probs[[s]]) <- paste0("pr_state_",1:m)
    }

    # # Sample states
    # state_seq[1:n_vary[s], s] <- apply(probs[[s]], 2, which.max)

    # for(h in 1:forecast_steps){
    #
    # }
    # colnames(probs[[s]]) <- paste0("pr_state_",1:m)
    # probs[[s]] <- cbind(subj = s, state = apply(probs[[s]], 1, which.max), horizon = 1:forecast_steps, probs[[s]])

    # Structure output
    colnames(probs[[s]]) <- paste0("pr_state_",1:m)
    if(return_all){
      probs[[s]] <- cbind(subj = s, state = apply(probs[[s]], 1, which.max), horizon = c(rep(0, n_vary[s]), 1:forecast_steps), probs[[s]])
    } else {
      probs[[s]] <- cbind(subj = s, state = apply(probs[[s]], 1, which.max), horizon = 1:forecast_steps, probs[[s]])
    }

    # # If a range of observed values is given, obtain the forecasting distribution (forecast probabilities):
    # if(!is.null(x)){
    #   for(q in 1:n_dep){
    #     # obs_probs[[q]][h,] <- emiss
    #
    #     # Forecast probs for scoring values
    #     probs[[1]][,4:6] %*% t(outer(seq(0,100,25), Y = emiss[[q]][,1], FUN = stats::dnorm, sd = emiss[[q]][,2]))
    #
    #     # Forecast prob of scoring over value
    #     probs[[1]][,4:6] %*% t(1-outer(seq(0,100,25), Y = emiss[[q]][,1], FUN = stats::pnorm, sd = emiss[[q]][,2]))
    #
    #     # Forecast prob of scoring under value
    #     probs[[1]][,4:6] %*% t(outer(seq(0,100,25), Y = emiss[[q]][,1], FUN = stats::pnorm, sd = emiss[[q]][,2]))
    #
    #   }
    #   # allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "categorical")
    #
    # }

  }

  forecast_probs <- do.call(rbind, lapply(probs, function(s) s))

  # state_seq <- as.data.frame(state_seq)
  # colnames(state_seq) <- id
  # state_seq_stacked <- utils::stack(state_seq)[,c(2,1)]
  # colnames(state_seq_stacked) <- c("subj", "state")
  # if(return_state_prob == TRUE){
  #   probs <- sapply(probs, t, simplify = FALSE)
  #   probs <- do.call(rbind, probs)
  #   colnames(probs) <- paste0("pr_state_", 1:m)
  #   state_seq_stacked <- cbind(state_seq_stacked, probs)
  # }
  # return(state_seq = state_seq_stacked)

  # Return results
  return(forecast_probs = forecast_probs)

}


#' Forecasts (forward predictions) using a multilevel hidden Markov model
#'
#' forecast_mHMM2 makes forecasts using a Bayesian one-step-look ahead
#' procedure, sampling states and observations iteratively over the forecasting
#' horizon range over J iterations of the MCMC sampler.
#'
#'
#' @export

# forecast_mHMM2 <- function(s_data, object, forecast_steps = 1, value_range = NULL, n_iter = NULL, burn_in = NULL, show_progress = TRUE) {
#
#   if (!is.mHMM(object)){
#     stop("The input object used should be from the class mHMM, obtained by using the function mHMM.")
#   }
#   if(sum(objects(object$PD_subj[[1]]) %in% "log_likl") != 1){
#     stop("The input object is created using an earlier version of the mHMMbayes package.
#          Please re-run the function mHMM with the current package version, or post-process the object
#          using the earlier version of the package.")
#   }
#   input      <- object$input
#   data_distr <- input$data_distr
#   id         <- unique(s_data[,1])
#   n_subj     <- length(id)
#   if(length(object$PD_subj) != n_subj){
#     stop("s_data used should be from the same subjects used for creating the object in mHMM.
#          The number of subjects in the datasets are not the same.")
#   }
#   n_vary     <- table(s_data[,1])
#   max_n      <- max(n_vary)
#   state_seq  <- matrix(NA, ncol = n_subj, nrow = max_n)
#   probs      <- vector(mode = "list", length = n_subj)
#   n_dep      <- input$n_dep
#   m          <- input$m
#   if(is.null(burn_in)){
#     burn_in  <- input$burn_in
#   }
#   J          <- n_iter <- input$J
#   if (burn_in >= (J-1)){
#     stop(paste("The specified burn in period should be at least 2 points smaller
#                compared to the number of iterations J, J =", J))
#   }
#
#   # Get subject-specific emissions
#   if(data_distr == "categorical"){
#     q_emiss    <- input$q_emiss
#     int_est_emiss <- rep(list(lapply(q_emiss-1, dif_matrix, rows = m)), n_subj)
#     est_emiss  <- rep(list(lapply(q_emiss, dif_matrix, rows = m)), n_subj)
#     for(s in 1:n_subj){
#       for(j in 1:n_dep){
#         int_est_emiss[[s]][[j]][] <- matrix(apply(object$emiss_int_subj[[s]][[j]][burn_in:J, ], 2, median),
#                                             byrow = TRUE, ncol = q_emiss[j]-1, nrow = m)
#         est_emiss[[s]][[j]][] <- int_to_prob(int_est_emiss[[s]][[j]])
#       }
#     }
#   } else if(data_distr == "continuous"){
#     est_emiss  <- rep(list(rep(list(matrix(NA_real_,nrow = m, ncol = 2)),n_dep)), n_subj)
#     for(s in 1:n_subj){
#       for(q in 1:n_dep){
#         est_emiss[[s]][[q]][] <- matrix(round(c(apply(object$PD_subj[[s]]$cont_emiss[((burn_in + 1): J),((q-1) * m + 1):(q * m)], 2, median),
#                                                 apply(object$PD_subj[[s]]$cont_emiss[((burn_in + 1): J), (n_dep * m + (q-1) * m + 1):(n_dep * m + q * m)], 2, median)),3),
#                                         ncol = 2, nrow = m)
#       }
#     }
#   } else if(data_distr == "count"){
#     est_emiss  <- rep(list(rep(list(matrix(NA_real_,nrow = m, ncol = 1)),n_dep)), n_subj)
#     for(s in 1:n_subj){
#       for(q in 1:n_dep){
#         est_emiss[[s]][[q]][] <- matrix(round(apply(object$PD_subj[[s]]$count_emiss[((burn_in + 1): J),((q-1) * m + 1):(q * m)], 2, median),3),
#                                         ncol = 1, nrow = m)
#       }
#     }
#   }
#
#   # Get subject-specific transitions
#   est_gamma <- obtain_gamma(object, level = "subject")
#
#   # Obtain the forward probabilities and make forecasts
#   gamma_iter <- vector("list", n_subj)
#   # trans <- rep(list(vector("list", m)), n_subj)
#   drawn_state <- rep(list(matrix(NA, nrow = forecast_steps, ncol = n_iter)), n_subj)
#   drawn_obs <- rep(list(rep(list(matrix(NA, nrow = forecast_steps, ncol = n_iter)), n_dep)), n_subj)
#
#   itime <- proc.time()[3]
#   if(show_progress == TRUE){
#     cat("Progress of the Bayesian mHMM algorithm:", "\n")
#     pb <- utils::txtProgressBar(min = 1, max = n_subj, style = 3)
#   }
#
#   for(s in 1:n_subj){
#
#     emiss   <- est_emiss[[s]]
#     gamma   <- est_gamma[[s]]
#
#     if(data_distr == "categorical"){
#       probs[[s]]    <- cat_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
#                                                         matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
#                                               m = m, emiss = emiss, gamma = gamma, n_dep = n_dep, delta=NULL)[[1]][,(n_vary[s]):(n_vary[s]+forecast_steps)]
#       # colnames(probs[[s]]) <- paste0("pr_state_",1:m)
#     } else if(data_distr == "continuous"){
#       probs[[s]]    <- cont_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
#                                                          matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
#                                                m = m, emiss = emiss, n_dep = n_dep, gamma = gamma)[[1]][,(n_vary[s]):(n_vary[s]+forecast_steps)]
#       # colnames(probs[[s]]) <- paste0("pr_state_",1:m)
#     } else if(data_distr == "count"){
#       probs[[s]]    <- count_mult_fw_r_to_cpp(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
#                                                           matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
#                                                 m = m, emiss = emiss, n_dep = n_dep, gamma = gamma)[[1]][,(n_vary[s]):(n_vary[s]+forecast_steps)]
#       # colnames(probs[[s]]) <- paste0("pr_state_",1:m)
#     }
#
#     # colnames(probs[[s]]) <- paste0("pr_state_",1:m)
#     # probs[[s]] <- cbind(subj = s, state = apply(probs[[s]], 1, which.max), horizon = 1:forecast_steps, probs[[s]])
#
#     # MCMC sampling:
#     # probs      <- vector(mode = "list", length = n_subj)
#
#     for (iter in 1:n_iter){
#
#       gamma_iter[[s]] <- matrix(object$PD_subj[[1]]$trans_prob[iter,], nrow = m, byrow = TRUE)
#       state0 	<- sample(1:m, 1, prob = c(probs[[s]][,1]))
#
#       for(h in 1:forecast_steps){
#
#         if(h == 1){
#
#           drawn_state[[s]][h, iter] 	<- sample(1:m, 1, prob = (gamma_iter[[s]][state0,]))
#
#         } else {
#
#           # Using the forward probabilites, sample the state sequence in a backward manner.
#           # In addition, saves state transitions in trans, and conditional observations within states in cond_y
#           drawn_state[[s]][h,iter] 	  <- sample(1:m, 1, prob = (gamma_iter[[s]][drawn_state[[s]][h-1, iter],]))
#
#         }
#
#         # Sample observations for the emission distribution
#         # for (i in 1:m){
#
#         i <- drawn_state[[s]][h,iter]
#
#         # Sample observations
#         for(q in 1:n_dep){
#           drawn_obs[[s]][[q]][h,iter] <- rnorm(n = 1,
#                                                mean = object$PD_subj[[s]]$cont_emiss[iter,(i+(q-1)*n_dep)],
#                                                sd = object$PD_subj[[s]]$cont_emiss[iter,(m*n_dep+i+(q-1)*n_dep)])
#         }
#         # }
#
#       }
#
#     }
#
#     # Save forward probabilities
#     probs[[s]] <- t(probs[[s]])
#     colnames(probs[[s]]) <- paste0("pr_state_",1:m)
#     probs[[s]] <- cbind(subj = s, state = apply(probs[[s]], 1, which.max)[-1], horizon = 1:forecast_steps, probs[[s]][-1,])
#
#     # Save predicted states sampled
#     # drawn_state[[s]][h, iter]
#
#     # Save predicted observations sampled
#     # drawn_obs[[s]][h,q]
#
#     if(show_progress == TRUE){
#       utils::setTxtProgressBar(pb, s)
#     }
#
#   }
#
#   if(show_progress == TRUE){
#     close(pb)
#   }
#
#   forecast_probs <- do.call(rbind, lapply(probs, function(s) s))
#
#   # state_seq <- as.data.frame(state_seq)
#   # colnames(state_seq) <- id
#   # state_seq_stacked <- utils::stack(state_seq)[,c(2,1)]
#   # colnames(state_seq_stacked) <- c("subj", "state")
#   # if(return_state_prob == TRUE){
#   #   probs <- sapply(probs, t, simplify = FALSE)
#   #   probs <- do.call(rbind, probs)
#   #   colnames(probs) <- paste0("pr_state_", 1:m)
#   #   state_seq_stacked <- cbind(state_seq_stacked, probs)
#   # }
#   # return(state_seq = state_seq_stacked)
#
#   # Return results
#   return(list(forecast_probs = forecast_probs, predicted_states = drawn_state, predicted_obs = drawn_obs))
#
# }
