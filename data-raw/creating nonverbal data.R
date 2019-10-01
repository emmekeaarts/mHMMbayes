setwd("/Users/emmeke/surfdrive/R package HMM/mHMMbayes/data-raw")

# creating nonverbal_cov (subject level predictors)
pred <- read.csv(file = "dat_beh_predictors.csv")
head(pred)
colnames(pred)[1] <- "id"
pred <- pred[1:10,]
std.CDI <- (pred$CDI_change - mean(pred$CDI_change)) / sd(pred$CDI_change)
std.SCA <- (pred$SCARED_change - mean(pred$SCARED_change)) / sd(pred$SCARED_change)
pred2 <- cbind(pred, round(std.CDI,2), round(std.SCA,2))
pred3 <- pred2[1:10,-c(1:5, 7,8)]
colnames(pred3)<- c("diagnosis", "std_CDI_change", " std_SCA_change")
pred3[,1] <- pred3[,1] - 1
nonverbal_cov <- pred3
library(devtools)
usethis::use_data(nonverbal_cov, overwrite = TRUE)

# creating nonverbal (patient and therapist nonverbal behavior over time), using first 10 couples
id <- pred[,1]
input1 <- paste("cluster data therapist ", id, ".csv", sep = "")
input2 <- paste("cluster data patient ", id, ".csv", sep = "")

HMM_dyad_data <- vector("list", 10)
workdata.pat <- numeric(1)
workdata.ther <- numeric(1)

for(i in 1:10){
  workdata.ther <- read.csv2(input1[i], header = TRUE, dec = ".", sep = ";")[,2:3]
  colnames(workdata.ther) <- paste("T", colnames(workdata.ther), sep = "_")
  # make sure factors levels are in correct order
  workdata.ther$T_Vocalizing <- factor(workdata.ther$T_Vocalizing, levels = c( "Voc_no", "Voc_speak", "Voc_bch"))
  workdata.ther$T_Vocalizing <- unclass(workdata.ther$T_Vocalizing)
  workdata.ther$T_Looking <- factor(workdata.ther$T_Looking, levels = c("Look_no","Look_look"))
  workdata.ther$T_Looking <- unclass(workdata.ther$T_Looking)
  # workdata.ther$T_Head <- factor(workdata.ther$T_Head, levels = c("Head_no","Head_nod", "Head_shake", "Head_other"))
  # workdata.ther$T_Head <- unclass(workdata.ther$T_Head)
  # workdata.ther$T_Leg <- factor(workdata.ther$T_Leg, levels = c( "Leg_no", "Leg_move"))
  # workdata.ther$T_Leg <- unclass(workdata.ther$T_Leg)
  # workdata.ther$T_Hand_L <- factor(workdata.ther$T_Hand_L, levels = c("Hand_L.no", "Hand_L.gest", "Hand_L.botou", "Hand_L.obtou", "Hand_L.other"))
  # workdata.ther$T_Hand_L <- unclass(workdata.ther$T_Hand_L)
  # workdata.ther$T_Hand_R <- factor(workdata.ther$T_Hand_R, levels = c("Hand_R.no", "Hand_R.gest", "Hand_R.botou", "Hand_R.obtou", "Hand_R.other"))
  # workdata.ther$T_Hand_R <- unclass(workdata.ther$T_Hand_R)

  workdata.pat <- read.csv2(input2[i], header = TRUE, dec = ".", sep = ";")[,2:3]
  colnames(workdata.pat) <- paste("P", colnames(workdata.pat), sep = "_")
  # make sure factors levels are in correct order
  workdata.pat$P_Vocalizing <- factor(workdata.pat$P_Vocalizing, levels = c("Voc_no", "Voc_speak", "Voc_bch"))
  workdata.pat$P_Vocalizing <- unclass(workdata.pat$P_Vocalizing)
  workdata.pat$P_Looking <- factor(workdata.pat$P_Looking, levels = c("Look_no", "Look_look"))
  workdata.pat$P_Looking <- unclass(workdata.pat$P_Looking)
  # workdata.pat$P_Head <- factor(workdata.pat$P_Head, levels = c("Head_no", "Head_nod", "Head_shake", "Head_other"))
  # workdata.pat$P_Head <- unclass(workdata.pat$P_Head)
  # workdata.pat$P_Leg <- factor(workdata.pat$P_Leg, levels = c("Leg_no", "Leg_move"))
  # workdata.pat$P_Leg <- unclass(workdata.pat$P_Leg)
  # workdata.pat$P_Hand_L <- factor(workdata.pat$P_Hand_L, levels = c("Hand_L.no", "Hand_L.gest", "Hand_L.botou", "Hand_L.obtou", "Hand_L.other"))
  # workdata.pat$P_Hand_L <- unclass(workdata.pat$P_Hand_L)
  # workdata.pat$P_Hand_R <- factor(workdata.pat$P_Hand_R, levels = c("Hand_R.no", "Hand_R.gest", "Hand_R.botou", "Hand_R.obtou", "Hand_R.other"))
  # workdata.pat$P_Hand_R <- unclass(workdata.pat$P_Hand_R)

  HMM_dyad_data[[i]] <- list(y = as.matrix(cbind(workdata.pat,workdata.ther), rownames.force = FALSE))

}
length(HMM_dyad_data)

nonverbal <- cbind(1,HMM_dyad_data[[1]]$y)
for(i in 2:10){
  nonverbal <- rbind(nonverbal,cbind(i,HMM_dyad_data[[i]]$y))
}
colnames(nonverbal) <- c("id", "p_vocalizing", "p_looking", "t_vocalizing", "t_looking")
library(devtools)
usethis::use_data(nonverbal, overwrite = TRUE)



