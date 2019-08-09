## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = FALSE----------------------------------------------------
library(mHMMbayes)
nonverbal <- data.frame(nonverbal)
head(nonverbal)

## ---- fig.width = 7.2, fig.height = 4, echo = FALSE----------------------
# set labels and colors for the observed behavioral categorical outcomes
Voc_lab <- c("Speaking", "Back channeling", "Not Speaking")
Look_lab <-  c("Looking", "Not looking")
Voc_col <- c("darkslategray3", "darkslategray4", "gray85")
Look_col <- c("goldenrod1", "gray85")
cols = list(Voc_col, Look_col, Voc_col, Look_col)

time_s  <- seq(1,900)
couple1 <- cbind(nonverbal[nonverbal$id == 1,], time_s)

par(mar = c(4.3, 6.6, 2.1, 1.1))
plot(x = 1, xlim = c(0,300), ylim = c(0.5,6), type = "n", las = 1, xlab = "Time in minutes", xaxt = "n", yaxt = "n", ylab = "")
axis(2, at = seq(1,4), tick = FALSE, labels = c("P_vocalizing", "P_Looking", "T_vocalizing", "T_Looking"), las = 1)
axis(1, at = seq(0,300,60), tick = TRUE, las = 1, labels = FALSE)
axis(1, at = seq(0,300,60), tick = FALSE, las = 1, labels = seq(1,6,1))
abline(v = seq(0,300,60), col = "gray85")

for(j in 2:5){
  for(i in 1:max(nonverbal[,j])){
    points(x = couple1$time_s[1:300][couple1[1:300,j] == i], 
           y = rep(j-1, sum(couple1[1:300,j] == i)), 
           pch = "|", col = cols[[j-1]][i])
  }
}

legend("topright", bty = "n", fill = Voc_col, legend = Voc_lab)
legend("topleft", bty = "n", fill = Look_col, legend = Look_lab)


## ---- include = FALSE----------------------------------------------------
# specifying general model properties:
m <- 2
n_dep <- 4
q_emiss <- c(3, 2, 3, 2)

# specifying starting values
start.TM <- diag(.8, m)
start.TM[lower.tri(start.TM) | upper.tri(start.TM)] <- .2
start.EM <- list(matrix(c(0.9, 0.05, 0.05, 
                          0.05, 0.05, 0.9), byrow = TRUE,
                         nrow = m, ncol = q_emiss[1]), # vocalizing patient
                  matrix(c(0.9, 0.1, 
                           0.9, 0.1), byrow = TRUE, nrow = m,
                         ncol = q_emiss[2]), # looking patient
                  matrix(c(0.05, 0.05, 0.9, 
                           0.9, 0.05, 0.05), byrow = TRUE,
                         nrow = m, ncol = q_emiss[3]), # vocalizing therapist
                  matrix(c(0.9, 0.1, 
                           0.9, 0.1), byrow = TRUE, nrow = m,
                         ncol = q_emiss[4])) # looking therapist

load("nonv_2st_1000it.rda")

## ---- eval = FALSE-------------------------------------------------------
#  library(mHMMbayes)
#  # specifying general model properties:
#  m <- 2
#  n_dep <- 4
#  q_emiss <- c(3, 2, 3, 2)
#  
#  # specifying starting values
#  start.TM <- diag(.8, m)
#  start.TM[lower.tri(start.TM) | upper.tri(start.TM)] <- .2
#  start.EM <- list(matrix(c(0.9, 0.05, 0.05, 0.05, 0.05, 0.9), byrow = TRUE,
#                           nrow = m, ncol = q_emiss[1]), # vocalizing patient
#                    matrix(c(0.9, 0.1, 0.9, 0.1), byrow = TRUE, nrow = m,
#                           ncol = q_emiss[2]), # looking patient
#                    matrix(c(0.05, 0.05, 0.9, 0.9, 0.05, 0.05), byrow = TRUE,
#                           nrow = m, ncol = q_emiss[3]), # vocalizing therapist
#                    matrix(c(0.9, 0.1, 0.9, 0.1), byrow = TRUE, nrow = m,
#                           ncol = q_emiss[4])) # looking therapist

## ---- eval = FALSE-------------------------------------------------------
#  # Run a model without covariate(s) and default priors:
#  set.seed(241632)
#  out1 <- mHMM_mnl(s_data = nonverbal, gen = list(m = m, n_dep = n_dep,
#                   q_emiss = q_emiss), start_val = list(start.TM, start.EM),
#                   mcmc = list(J = 1000, burn_in = 200))
#  

## ------------------------------------------------------------------------
out1

## ------------------------------------------------------------------------
summary(out1)

## ------------------------------------------------------------------------
# When not specified, level defaults to "group"
gamma_pop <- obtain_gamma(out1)
gamma_pop

# To obtain the subject specific parameter estimates:
gamma_subj <- obtain_gamma(out1, level = "subject")
gamma_subj

## ---- fig.show='hold'----------------------------------------------------
# Transition probabilities at the group level and for subject number 1, respectively:
plot(gamma_pop)
plot(gamma_subj, subj_nr = 1)

## ---- fig.width = 7.2, fig.height = 4------------------------------------
# required input information for the plot
Voc_col <- c("darkslategray3", "darkslategray4", "gray85")
Voc_lab <- c("Speaking", "Back-channeling", "Not speaking")
burn_in <- 200
J <- 1000
qp <- 1				# dependent variable we are looking at
start <- c(0, q_emiss * m) 
start2 <- c(0, seq(from = (q_emiss[qp]-1) * 2, to = (q_emiss[qp]-1) * 2 * m, by = (q_emiss[qp]-1) * 2))

par(mfrow = c(1,2))
for(i in 1:m){
  # determining the scale of the y axis  
  max <- 0
    for(q in 1:q_emiss[qp]){
      new <- max(density(out1$emiss_prob_bar[[qp]][burn_in:J,q_emiss[qp] * (i-1) + q])$y)
      if(new > max){max <- new}	
    }
  # set plotting area  
  plot(x = 1, ylim = c(0, max), xlim = c(0,1), type = "n", 
       main = paste("Patient vocalizing, state", i), 
       yaxt = "n", ylab = "Density", xlab = "Conditional probability")
  for(q in 1:q_emiss[qp]){
    # add density curve for population level posterior distribution
    lines(density(out1$emiss_prob_bar[[qp]][burn_in:J,q_emiss[qp] * (i-1) + q]), 
          type = "l", col = Voc_col[q], lwd = 2.5)
    # add density curves for subject posterior distributions
    for(s in 1:10){
      lines(density(out1$PD_subj[[s]][burn_in:J,(sum(start[1:qp]) 
                                                 + (i-1)*q_emiss[qp] + q)]), 
            type = "l", col = Voc_col[q], lwd = 1.5, lty = 3)
    }	
  }
  legend("topright", col = Voc_col, legend = Voc_lab, bty = 'n', lty = 1, lwd = 2, cex = .7)
}

## ---- include = FALSE----------------------------------------------------
load("nonv_3st_1000it.rda")
load("nonv_4st_1000it.rda")

## ---- fig.width = 5, fig.height = 3--------------------------------------
summary(out3)
plot(obtain_gamma(out3))

