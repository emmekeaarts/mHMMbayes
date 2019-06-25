


plot_gamma <- function(object, level = "group", subj = 1){
  if (!requireNamespace("alluvial", quietly = TRUE)) {
    stop("Package \"alluvial\" needed for this function to work. Please install it.",
         call. = FALSE)
  }


  alluvial::alluvial()

}


# vignette("alluvial")

# sankey diagram (also known as alluvial diagram, or riveerplots. We use the R package alluvial)

From <- paste("State", rep(c(1,2,3), each = 3))
To <-  paste("State", rep(c(1,2,3), 3))

trans <- matrix(c(.6, .3, .1,
                  .2, .8, .0,
                  .3, .3, .4), byrow = TRUE, ncol = 3)
trans2 <- as.vector(t(trans))

out <- data.frame(From, To, trans2)

alluvial(out[,1:2], freq=out$trans2,
         cex = .8,
         col = c(rep(c("green","red", "blue"), eac = 3)),
         hide = out$trans2 == 0)
