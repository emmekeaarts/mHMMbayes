
plot.mHMM <- function(x, component = "gamma", dep = 1, col, cat_lab,
                      dep_lab, lwd1 = 2, lwd2 = 1, lty1 = 1, lty2 = 3,
                      legend_cex, burn_in, ...){
  if (!is.mHMM(x)){
    stop("The input object x should be from the class mHMM, obtained with the function mHMM.")
  }
  if (component != "gamma" & component != "emiss"){
    stop("The input specified under component should be a string, restrectid to state either gamma or emiss.")
  }
  object <- x
  input   <- x$input
  n_subj  <- input$n_subj
  if (missing(burn_in)){
    burn_in <- input$burn_in
  }
  J       <- input$J
  if (burn_in >= (J-1)){
    stop(paste("The specified burn in period should be at least 2 points smaller
               compared to the number of iterations J, J =", J))
  }
  old_par <- graphics::par(no.readonly =TRUE)
  on.exit(graphics::par(old_par))
  m       <- input$m
  q_emiss <- input$q_emiss
  n_dep   <- input$n_dep

  # if plotting gamma
  if(component == "gamma"){
    # if ggplot2 and dplyr is available
    if(nzchar(system.file(package = "ggplot2")) && nzchar(system.file(package = "dplyr"))){
      # group level gamma
      gg <- object$gamma_prob_bar %>% tibble::as_tibble() %>% dplyr::slice(burn_in:J) %>% dplyr::mutate(id = "group")
      # extract the peak of density for group level gamma
      peak <- gg %>% dplyr::select(where(is.numeric)) %>% apply(2, function(x) stats::density(x)$y) %>% max()
      # subject level gamma
      sg <- object$PD_subj %>% purrr::map(
        ~ tibble::as_tibble(.x) %>% dplyr::select(tidyselect::starts_with("S")) %>% dplyr::slice(burn_in:J)
      ) %>% dplyr::bind_rows(.id="id")
      # structure df for plotting
      df <- rbind(sg, gg) %>%
        tidyr::pivot_longer(cols = !id, names_to = "class", values_to = "value") %>%
        mutate(from_state = stringr::str_split(class, "to", simplify = T)[,1],
               to_state = stringr::str_split(class, "to", simplify = T)[,2],
               lty = as.factor(ifelse(id == "group", 1, 2)))
      # create plots
      plt <- df %>%
        ggplot2::ggplot(ggplot2::aes(x = value, color = to_state, alpha = id, linetype = lty, linewidth = lty)) +
        ggplot2::ylim(0, peak) +
        ggplot2::geom_density() +
        ggplot2::scale_linewidth_manual(values = c("1" = 0.7, "2" = 0.3))+
        ggplot2::scale_linetype_discrete(name = "level", labels = c("group", "subject"))+
        ggplot2::scale_color_discrete(name = "", labels = paste0("to S", 1:m))+
        ggplot2::facet_grid(~from_state, labeller = as_labeller(function(string, prefix = "from") paste(prefix, string))) +
        ggplot2::guides(alpha = "none", size = "none") +
        ggplot2::theme_bw()+
        ggplot2::theme(panel.spacing.x = unit(4, "mm"),
                       legend.position = "bottom",
                       axis.text.y = element_blank()
        ) +
        ggplot2::labs(x = "transition probability")
      ggplot2:::print.ggplot(plt) %>% suppressWarnings()

      # if ggplot2 is not available, original function follows
    } else {
      if (missing(col)){
        state_col <- grDevices::rainbow(m)
      } else {
        state_col <- col
      }
      if(m > 3){
        graphics::par(mfrow = c(2,ceiling(m/2)), mar = c(4,2,3,1) + 0.1, mgp = c(2,1,0))
      } else {
        graphics::par(mfrow = c(1,m), mar = c(4,2,3,1) + 0.1, mgp = c(2,1,0))
      }
      for(i in 1:m){
        max <- 0
        for(j in 1:m){
          new <- max(stats::density(object$gamma_prob_bar[burn_in:J, m * (i-1) + j])$y)
          if(new > max){max <- new}
        }
        graphics::plot.default(x = 1, ylim = c(0, max), xlim = c(0,1), type = "n", cex = .8,  main =
                                 paste("From state", i, "to state ..."), yaxt = "n", ylab = "",
                               xlab = "Transition probability", ...)
        graphics::title(ylab="Density", line=.5)
        for(j in 1:m){
          graphics::lines(stats::density(object$gamma_prob_bar[burn_in:J,m * (i-1) + j]),
                          type = "l", col = state_col[j], lwd = lwd1, lty = lty1)
          for(s in 1:n_subj){
            graphics::lines(stats::density(object$PD_subj[[s]][burn_in:J,(sum(q_emiss * m) + m * (i-1) + j)]),
                            type = "l", col = state_col[j], lwd = lwd2, lty = lty2)
          }
        }
        graphics::legend("topright", col = state_col, legend = paste("To state", 1:m),
                         bty = 'n', lty = 1, lwd = 2, cex = .8)
      }
    }

    # if plotting emiss
  } else if (component == "emiss"){
    if (missing(cat_lab)){
      cat_lab <- paste("Category", 1:q_emiss[dep])
    }
    if (missing(dep_lab)){
      dep_lab <- input$dep_labels[dep]
    }
    # if ggplot2 and dplyr is available
    if(nzchar(system.file(package = "ggplot2")) && nzchar(system.file(package = "dplyr"))){
      # group level emission
      ge <- object$emiss_prob_bar[[dep]] %>% tibble::as_tibble() %>% dplyr::slice(burn_in:J) %>% dplyr::mutate(id = "group") %>% dplyr::relocate(id)
      # extract the peak
      peak <- ge %>% dplyr::select(where(is.numeric)) %>% apply(2, function(x) stats::density(x)$y) %>% max()
      # subject level emission
      se <- object$PD_subj %>% purrr::map(
        ~ tibble::as_tibble(.x) %>% dplyr::slice(burn_in:J) %>% dplyr::select(tidyselect::contains(paste0("q", dep, collapse="")))
      )  %>% dplyr::bind_rows(.id="id") %>% magrittr::set_colnames(colnames(ge))

      # structure df for plotting
      df <- rbind(se, gg) %>%
        tidyr::pivot_longer(cols = !id, names_to = "class", values_to = "value") %>%
        dplyr::mutate(category = stringr::str_split(class, "_", simplify = T)[,1],
               state = stringr::str_split(class, "_", simplify = T)[,2],
               lty = as.factor(ifelse(id == "group", 1, 2)))

      plt <- df %>%
        ggplot2::ggplot(ggplot2::aes(x = value, color = category, alpha = id, linetype = lty, linewidth = lty)) +
        ggplot2::geom_density() +
        ggplot2::ylim(0, peak) +
        ggplot2::scale_linewidth_manual(values = c("1" = 0.5, "2" = 0.3))+
        ggplot2::scale_linetype_discrete(name = "level", labels = c("group", "subject"))+
        ggplot2::scale_color_discrete(name = "", labels = cat_lab) +
        ggplot2::facet_grid(~state, labeller = as_labeller(function(string, prefix = paste0(dep_lab, ", ")) paste(prefix, string))) +
        ggplot2::guides(alpha = "none", size = "none") +
        ggplot2::theme_bw()+
        ggplot2::theme(panel.spacing.x = unit(4, "mm"),
              legend.position = "bottom",
              axis.text.y = element_blank()
        ) +
        ggplot2::labs(x = "conditional probability")

      ggplot2:::print.ggplot(plt) %>% suppressWarnings()

      # if ggplot2 is not available, original function follows
    } else {
      start <- c(0, q_emiss * m)
      start2 <- c(0, seq(from = (q_emiss[dep]-1) * 2, to = (q_emiss[dep]-1) * 2 * m, by = (q_emiss[dep]-1) * 2))
      if (missing(col)){
        cat_col <- grDevices::rainbow(q_emiss[dep])
      } else {
        cat_col <- col
      }
      if(m > 3){
        graphics::par(mfrow = c(2,ceiling(m/2)), mar = c(4,2,3,1) + 0.1, mgp = c(2,1,0))
      } else {
        graphics::par(mfrow = c(1,m), mar = c(4,2,3,1) + 0.1, mgp = c(2,1,0))
      }
      for(i in 1:m){
        # determining the scale of the y axis
        max <- 0
        for(q in 1:q_emiss[dep]){
          new <- max(stats::density(object$emiss_prob_bar[[dep]][burn_in:J,q_emiss[dep] * (i-1) + q])$y)
          if(new > max){max <- new}
        }
        # set plotting area
        graphics::plot.default(x = 1, ylim = c(0, max), xlim = c(0,1), type = "n",
                               main = paste(dep_lab, ", state", i),
                               yaxt = "n", ylab = "", xlab = "Conditional probability", ...)
        graphics::title(ylab="Density", line=.5)
        for(q in 1:q_emiss[dep]){
          # add density curve for population level posterior distribution
          graphics::lines(stats::density(object$emiss_prob_bar[[dep]][burn_in:J,q_emiss[dep] * (i-1) + q]),
                          type = "l", col = cat_col[q], lwd = lwd1, lty = lty1)
          # add density curves for subject posterior distributions
          for(s in 1:n_subj){
            graphics::lines(stats::density(object$PD_subj[[s]][burn_in:J,(sum(start[1:dep])
                                                                          + (i-1)*q_emiss[dep] + q)]),
                            type = "l", col = cat_col[q], lwd = lwd2, lty = lty2)
          }
        }
        graphics::legend("topright", col = cat_col, legend = cat_lab, bty = 'n', lty = 1, lwd = 2, cex = .7)
      }
    }

  }
}






