

plot.mHMM_gamma <- function(x, subj_nr = NULL, cex = 0.8, col, hide, ...){
  if (!requireNamespace("alluvial", quietly = TRUE)) {
    stop("Package \"alluvial\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Package \"grDevices\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!is.mHMM_gamma(x)){
    stop("The input object x should be from the class mHMM_gamma, obtained with the function obtain_gamma.")
  }
  old_par <- graphics::par(no.readonly =TRUE)
  on.exit(graphics::par(old_par))
  if (is.list(x)){
    if (is.null(subj_nr)){
      stop("When the input object x represents the subject specific transition
           probability matrices, the subject for which the probabilities should
           be plotted needs to be specified with the input variable -subj_nr-.")
    }
    m <- dim(x[[subj_nr]])[1]
    From <- paste("State", rep(1:m, each = m))
    To <-  paste("State", rep(1:m, m))
    trans <- as.vector(t(x[[subj_nr]]))
    foo <- data.frame(From, To, trans)
    if(missing(col)){
      col <- c(rep(grDevices::rainbow(m), each = m))
    }
    if (missing(hide)){
      hide <- foo$trans == 0
    }
    # if ggplot2 and ggalluvial is available
    if(nzchar(system.file(package = "ggplot2")) && nzchar(system.file(package = "ggalluvial"))){
      ggplot2::ggplot(foo, ggplot2::aes(axis1 = From, axis2 = To, y = trans)) +
        ggalluvial::geom_alluvium(ggplot2::aes(fill = From)) +
        ggalluvial::geom_stratum() +
        ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = ggplot2::after_stat(stratum))) +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "bottom",
                       plot.title = element_text(hjust = 0.5)) +
        ggplot2::labs(title= paste0("Transition probabilities for subject", subj_nr))
      #if ggplot2 is not available, original function follows
    } else {
    alluvial::alluvial(foo[,1:2], freq=foo$trans,
                       cex = cex,
                       col = col,
                       hide = hide, ...)
    }
  } else {
    if(!is.null(subj_nr)){
      warning("The subject number can only be specified when plotting the subject level transition probabilities. Currently, the group level transition probabilities are plotted.")
    }
    m <- dim(x)[1]
    From <- paste("State", rep(1:m, each = m))
    To <-  paste("State", rep(1:m, m))
    trans <- as.vector(t(x))
    foo <- data.frame(From, To, trans)
    if(missing(col)){
      col <- c(rep(grDevices::rainbow(m), each = m))
    }
    if (missing(hide)){
      hide <- foo$trans == 0
    }
    # if ggplot2 and ggalluvial is available
    if(nzchar(system.file(package = "ggplot2")) && nzchar(system.file(package = "ggalluvial"))){
      ggplot2::ggplot(foo, ggplot2::aes(axis1 = From, axis2 = To, y = trans)) +
        ggalluvial::geom_alluvium(ggplot2::aes(fill = From)) +
        ggalluvial::geom_stratum() +
        ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = ggplot2::after_stat(stratum))) +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "bottom",
                       plot.title = element_text(hjust = 0.5)) +
        ggplot2::labs(title= "Transition probabilities at the group level")
      #if ggplot2 is not available, original function follows
    } else {
    alluvial::alluvial(foo[,1:2], freq=foo$trans,
                       cex = cex,
                       col = col,
                       hide =  hide, ...)
    }
  }
}


