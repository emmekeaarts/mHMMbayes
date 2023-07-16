

# w = xx_emiss[[1]][,2]
# covariate <- w

v = sample(c(0,1), object$input$n_subj, replace = T)
covariate <- v

#' @param covariate A numeric vector specifying the values of a single covariate.
plot_pred <- function(object, component = "gamma", covariate, dep = 1, cat_lab, dep_lab, ...) {
  input   <- object$input
  dep_labels <- input$dep_labels
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  q_emiss <- input$q_emiss
  n_dep   <- input$n_dep
  covar_type_gamma <- input$covar_type[[1]]
  covar_type_emiss <- input$covar_type[[2]] # given that all DVs share the same covariates

  # if plotting gamma
  if(component == "gamma"){
    # prepare subject df
    gg_subj <- object$PD_subj %>%
      purrr::imap(\(x, idx)  tibble::as_tibble(x) %>%
                    dplyr::mutate(covariate = covariate[idx]) %>%
                    dplyr::select(tidyselect::starts_with("S"), covariate) %>%
                    dplyr::slice(burn_in:J)
      ) %>%
      dplyr::bind_rows(.id="subject") %>%
      tidyr::pivot_longer(cols=!c(subject, covariate), names_to = "states", values_to="prob") %>%
      dplyr::mutate(From_state = stringr::str_split(states, "to", simplify = T)[,1],
                    To_state = stringr::str_split(states, "to", simplify = T)[,2],
                    dplyr::across(From_state:To_state, stringr::str_replace, "S", "State ")
      )

    # prepare emiss predicted value df
    preds <- pred_probs(object, covariate = covariate) %>%
      tidyr::pivot_longer(!c(From_state, covariate), names_to = "To_state", values_to = "prob") %>%
      dplyr::mutate(dplyr::across(To_state, stringr::str_replace, "To_state_", "State "))

    plot <- if(covar_type_gamma == "continuous") {
      ggplot2::ggplot(preds, ggplot2::aes(x = covariate, y = prob, color = To_state)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = gg_subj, ggplot2::aes(x = covariate, y = prob, color = To_state), alpha = 0.1, size = 0.3) +
      ggplot2::labs(x = "Covariate", y = "Transition probability", title = "Predicted transition probability per covariate value", color = "To state") +
      ggplot2::ylim(0,1) +
      ggplot2::facet_grid(~From_state, labeller = as_labeller(function(string, prefix = "From") paste(prefix, string))) +
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = "bottom")

    } else if(covar_type_gamma == "dichotomous") {
      ggplot2::ggplot(preds, ggplot2::aes(x = factor(covariate), y = prob, color = To_state)) +
        ggplot2::geom_boxplot(outlier.size = 0.5, outlier.alpha=0.3, size = 0.3) +
        ggplot2::geom_jitter(data = gg_subj, ggplot2::aes(x = factor(covariate), y = prob, color = To_state), alpha = 0.1, size = 0.3) +
        ggplot2::labs(x = "Covariate", y = "Transition probability", title = "Predicted transition probability per covariate value", color = "To state") +
        ggplot2::ylim(0,1) +
        ggplot2::facet_grid(~From_state, labeller = as_labeller(function(string, prefix = "From") paste(prefix, string))) +
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position = "bottom")

    }

  }

  # if plotting emiss
  else if(component == "emiss"){
    if (missing(cat_lab)){
      cat_lab <- paste("Category", 1:q_emiss[dep])
    }
    if (missing(dep_lab)){
      dep_lab <- input$dep_labels[dep]
    }
    # prepare subject df
    gg_subj <- object$PD_subj %>%
      purrr::imap(\(x, idx)  tibble::as_tibble(x) %>%
                    dplyr::mutate(covariate = covariate[idx]) %>%
                    dplyr::select(tidyselect::contains(paste0("q", dep, collapse="")), covariate) %>%
                    dplyr::slice(burn_in:J)
      ) %>%
      dplyr::bind_rows(.id = "subject") %>%
      tidyr::pivot_longer(cols=!c(subject, covariate), names_to = "cats", values_to="prob") %>%
      dplyr::mutate(Category = stringr::str_split(cats, "_", simplify = T)[,2],
                    State = stringr::str_split(cats, "_", simplify = T)[,3],
                    dplyr::across(Category, stringr::str_replace, "emiss", "Category"),
                    dplyr::across(State, stringr::str_replace, "S", "State")
      )

    # prepare emiss predicted values
    preds <- pred_probs(object, covariate = covariate, component = "emiss") %>%
      tidyr::pivot_longer(!c(State, covariate), names_to = "Category", values_to = "prob")

    plot <- if(covar_type_emiss == "continuous"){
      ggplot2::ggplot(preds, ggplot2::aes(x = covariate, y = prob, color = Category)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = gg_subj, aes(x = covariate, y = prob, color = Category),
                 alpha = 0.1, size = 0.3) +
      ggplot2::labs(x = "Covariate", y = "Emission probability", title = "Predicted emission probability per covariate value") +
      ggplot2::ylim(0,1) +
      ggplot2::facet_grid(~State) +
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = "bottom")

    } else if(covar_type_emiss == "dichotomous") {
      ggplot2::ggplot(preds, ggplot2::aes(x = factor(covariate), y = prob, color = Category)) +
        ggplot2::geom_boxplot(outlier.size = 0.5, outlier.alpha=0.3, size = 0.3) +
        ggplot2::geom_jitter(data = gg_subj, aes(x = factor(covariate), y = prob, color = Category), alpha = 0.1, size = 0.3) +
        ggplot2::labs(x = "Covariate", y = "Emission probability", title = "Predicted emission probability per covariate value", color = "To state") +
        ggplot2::ylim(0,1) +
        ggplot2::facet_grid(~State) +
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position = "bottom")
    }
  }
  return(plot)
}


# plot_pred(object, covariate = w)
# plot_pred(object, covariate = v)
#
# plot_pred(object, covariate = w, component = "emiss")
# plot_pred(object, covariate = v, component = "emiss")
