

# w = xx_emiss[[1]][,2]
covariate <- w

v = sample(c(0,1), object$input$n_subj, replace = T)
covariate <- v

#' @param covariate A numeric vector specifying the values of a single covariate for all subjects (i.e., the length of vector should be equal to the number of subject).
#'
#' @export

plot_pred <- function(object, component = "gamma", covariate, dep = 1, cat_lab, dep_lab, ...) {
  input   <- object$input
  n_subj  <- input$n_subj
  dep_labels <- input$dep_labels
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  q_emiss <- input$q_emiss
  n_dep   <- input$n_dep
  covar_type_gamma <- input$covar_type[[1]]
  covar_type_emiss <- input$covar_type[[2]] # given that all DVs share the same covariates

  # common theme
  common_theme <- ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  # if plotting gamma
  if(component == "gamma"){
    # prepare subject df
    gg_subj <- obtain_gamma(object, level = "subject", burn_in = burn_in) |>
    purrr::reduce(rbind) |>
    tibble::as_tibble(rownames="From_state") |>
    dplyr::mutate(subject = rep(1:n_subj, each = m),
           covariate = rep(covariate, each = m)) |>
      tidyr::pivot_longer(cols=!c(subject, covariate, From_state),
                          names_to = "To_state", values_to="prob") |>
      dplyr::mutate(dplyr::across(c(To_state, From_state), stringr::str_replace,
                                  "To state|From state", "State"))

    # gg_subj <- object$gamma_int_subj |>
    #   purrr::imap(\(x, idx)  tibble::as_tibble(x) |>
    #                 dplyr::slice(-(1:burn_in)) |>
    #                 dplyr::summarize(across(everything(), median)) |>
    #                 dplyr::mutate(covariate = covariate[idx])
    #               )|>
    #   dplyr::bind_rows(.id="subject") |>
    #   tidyr::pivot_longer(cols=!c(subject, covariate), names_to = "states", values_to="prob") |>
    #   dplyr::mutate(From_state = paste("State", str_extract_all(states, "(\\d)", simplify=T)[,1]),
    #                 To_state = paste("State", str_extract_all(states, "(\\d)", simplify=T)[,2]))

    # prepare emiss predicted value df
    preds <- pred_probs(object, covariate = covariate) |>
      tidyr::pivot_longer(!c(From_state, covariate), names_to = "To_state", values_to = "prob") |>
      dplyr::mutate(dplyr::across(To_state, stringr::str_replace, "To_state_", "State "))

    if(covar_type_gamma == "continuous") {
      plot <- ggplot2::ggplot(data = gg_subj, ggplot2::aes(x = covariate, y = prob, color = To_state)) +
        ggplot2::geom_point(alpha = 0.6, size = 0.7) +
        ggplot2::geom_line(data = preds) +
        ggplot2::scale_color_brewer(palette="Accent") +
        ggplot2::labs(x = "Covariate", y = "Transition probability", title = "Predicted transition probability per covariate value", color = "To state") +
        ggplot2::ylim(0,1) +
        ggplot2::facet_grid(~From_state, labeller = ggplot2::as_labeller(function(string, prefix = "From") paste(prefix, string))) +
        common_theme

    } else if(covar_type_gamma == "dichotomous") {
      plot <- ggplot2::ggplot(data = gg_subj, ggplot2::aes(x = factor(covariate), y = prob, color = To_state)) +
        geom_point(alpha = 0.6, size = 0.7) +
        ggplot2::geom_boxplot(data = preds, width = 0.7, size = 0.5) +
        ggplot2::scale_color_brewer(palette="Accent") +
        ggplot2::labs(x = "Covariate", y = "Transition probability", title = "Predicted transition probability per covariate value", color = "To state") +
        ggplot2::ylim(0,1) +
        ggplot2::facet_grid(~From_state, labeller = ggplot2::as_labeller(function(string, prefix = "From") paste(prefix, string))) +
        common_theme
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
    gg_subj <- obtain_emiss(object, level = "subject", burn_in = burn_in)[[dep]] |>
      purrr::reduce(rbind) |>
      tibble::as_tibble(rownames="State") |>
      dplyr::mutate(subject = rep(1:n_subj, each = m),
                    covariate = rep(covariate, each = m)) |>
      tidyr::pivot_longer(cols=!c(subject, covariate, State),
                          names_to = "Category", values_to="prob")

    # gg_subj <- object$emiss_int_subj |>
    #   purrr::imap(\(x, idx)  tibble::as_tibble(x) |>
    #                 dplyr::mutate(covariate = covariate[idx]) |>
    #                 dplyr::select(tidyselect::contains(paste0("q", dep, collapse="")), covariate) |>
    #                 dplyr::slice(-(1:burn_in))
    #   ) |>
    #   dplyr::bind_rows(.id = "subject") |>
    #   tidyr::pivot_longer(cols=!c(subject, covariate), names_to = "cats", values_to="prob") |>
    #   dplyr::mutate(Category = stringr::str_split(cats, "_", simplify = T)[,2],
    #                 State = stringr::str_split(cats, "_", simplify = T)[,3],
    #                 dplyr::across(Category, stringr::str_replace, "emiss", "Category"),
    #                 dplyr::across(State, stringr::str_replace, "S", "State")
    #   )

    # prepare emiss predicted values
    preds <- pred_probs(object, covariate = covariate, component = "emiss", dep = dep) |>
      tidyr::pivot_longer(!c(State, covariate), names_to = "Category", values_to = "prob")

    if(covar_type_emiss == "continuous"){
      plot <- ggplot2::ggplot(data = gg_subj, ggplot2::aes(x = covariate, y = prob, color = Category)) +
        ggplot2::geom_point(alpha = 0.6, size = 0.7) +
        ggplot2::geom_line(data = preds) +
      ggplot2::scale_color_brewer(palette="Accent") +
      ggplot2::labs(x = "Covariate", y = "Emission probability", title = "Predicted emission probability per covariate value") +
      ggplot2::ylim(0,1) +
      ggplot2::facet_grid(~State) +
      common_theme

    } else if(covar_type_emiss == "dichotomous") {
      plot <- ggplot2::ggplot(data = gg_subj, ggplot2::aes(x = factor(covariate), y = prob, color = Category)) +
        ggplot2::geom_point(alpha = 0.6, size = 0.7) +
        ggplot2::geom_boxplot(data = preds, width = 0.7, size = 0.5) +
        ggplot2::scale_color_brewer(palette="Accent") +
        ggplot2::labs(x = "Covariate", y = "Emission probability", title = "Predicted emission probability per covariate value", color = "To state") +
        ggplot2::ylim(0,1) +
        ggplot2::facet_grid(~State) +
        common_theme
    }
  }
  return(plot)
}

utils::globalVariables(c("subject", "states", "prob", "cats", "Category"))


plot_pred(object, covariate = w)
plot_pred(object, covariate = v)

plot_pred(object, covariate = w, component = "emiss")
plot_pred(object, covariate = v, component = "emiss")
