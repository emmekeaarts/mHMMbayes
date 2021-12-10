library(mHMMbayes)
?mHMMbayes
vignette("estimation-mhmm")
data("nonverbal")
summary(nonverbal)
m <- 3 #number of states
n_dep <- 4 #number of dependent variables used to infer the hidden states
q_emiss <- c(3, 2, 3, 2) #number of categorical outcomes of each of dependent variables
start_TM <- diag(.8, m) #transition matrix
start_TM
start_TM[lower.tri(start_TM) | upper.tri(start_TM)] <- .2
start_EM <- list(matrix(c(0.05, 0.90, 0.05,
                          0.90, 0.05, 0.05,
                          0.60, 0.30, 0.10), byrow = TRUE,
                        nrow = m, ncol = q_emiss[1]), # vocalizing patient
                 matrix(c(0.1, 0.9,
                          0.1, 0.9,
                          0.5, 0.5), byrow = TRUE, nrow = m,
                        ncol = q_emiss[2]), # looking patient
                 matrix(c(0.90, 0.05, 0.05,
                          0.05, 0.90, 0.05,
                          0.60, 0.30, 0.10), byrow = TRUE,
                        nrow = m, ncol = q_emiss[3]), # vocalizing therapist
                 matrix(c(0.1, 0.9,
                          0.1, 0.9,
                          0.5, 0.5), byrow = TRUE, nrow = m,
                        ncol = q_emiss[4])) # looking therapist


set.seed(14532)
out_3st <- mHMM(s_data = nonverbal, gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),start_val = c(list(start_TM), start_EM),mcmc = list(J = 100, burn_in = 5))

group_emiss3<-obtain_emiss(out_3st,level = "group")

subject_emiss3<-obtain_emiss(out_2st,level = "subject")
data("nonverbal")
m1 <- 2 #number of states
n_dep1 <- 4 #number of dependent variables used to infer the hidden states
q_emiss1 <- c(3, 2, 3, 2) #number of categorical outcomes of each of dependent variables
start_TM1 <- diag(.8, m1) #transition matrix
start_TM1
start_TM1[lower.tri(start_TM) | upper.tri(start_TM)] <- .2
start_EM1 <- list(matrix(c(0.05, 0.90, 0.05,
                           0.90, 0.05, 0.05), byrow = TRUE,
                         nrow = m1, ncol = q_emiss1[1]), # vocalizing patient
                  matrix(c(0.1, 0.9,
                           0.1, 0.9), byrow = TRUE, nrow = m1,
                         ncol = q_emiss1[2]), # looking patient
                  matrix(c(0.90, 0.05, 0.05,
                           0.05, 0.90, 0.05), byrow = TRUE,
                         nrow = m1, ncol = q_emiss1[3]), # vocalizing therapist
                  matrix(c(0.1, 0.9,
                           0.1, 0.9), byrow = TRUE, nrow = m1,
                         ncol = q_emiss1[4])) # looking therapist


set.seed(14532)
out_2st <- mHMM(s_data = nonverbal, gen = list(m = m1, n_dep = n_dep1, q_emiss = q_emiss1),start_val = c(list(start_TM1), start_EM1),mcmc = list(J = 100, burn_in = 10))

group1<-obtain_emiss(out_2st,level = "group")

library(RColorBrewer)
coul <- brewer.pal(max(q_emiss), "Pastel2")


group<-c("voc_patient","look_patient","voc_therapist","look_therapist")

subject_emiss<-obtain_emiss(out_3st,level = "subject")
# advice on finding colours. Look
library(RColorBrewer)
n_col <- 60 #here you can add a number of cathegories you have got
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colr<-sample(col_vector, n_col) # and here it is an array of colors created



barplot_emiss<- function(emiss,q_emiss,m,n_dep,group_labs,color,nr_sub=NULL){
  #emiss is obtain_emiss object
  #q_emiss is vector of numbers of cathegories for each observed variable
  #m is number of states
  #n_dep number of dependent variables
  #group_labs are observation group names
  #color vector of colors and it has to be equal to the maximal number of cathegories
  #nr_sub if provided it is number of

  #number of colours has to be equal max(q_emiss)
  new<-list()
  if(is.null(nr_sub)==F){
    for(n in 1:n_dep){
      new[[n]]<-emiss[[n]][[nr_sub]]
    }
    emiss=new
  }

  emission_states<-list()
  data<-matrix(nrow=n_dep,ncol=max(q_emiss))
  for(s in 1:m){

    for(i in 1:n_dep){
      while(ncol(emiss[[i]])!=max(q_emiss)) emiss[[i]]<-cbind(emiss[[i]],0)
      data[i,]<-emiss[[i]][s,]
    }
    data<-as.data.frame(data)
    row.names(data)<-c(group_labs)
    emission_states[[s]]<-as.data.frame(data)
    data<-matrix(nrow=n_dep,ncol=max(q_emiss))
  }
  a<-emission_states
  par(mfrow=c(1,m))
  for(j in 1:m){
    if(j==1){
      if(is.null(nr_sub)==F){
        main_sub<-paste("Subject ",nr_sub)
        par(mar = c(6, 1, 4, 2))
        barplot(t(a[[j]]),xlab=paste("state",j),yaxt="n",legend=F, col=coul,las=2)
        title(main=paste(main_sub," emission distribution"), adj=0.1)
      }
      else{
        par(mar = c(6, 1, 4, 2))
        barplot(t(a[[j]]),main=paste("Group emission distribution"),xlab=paste("state",j),yaxt="n",legend=F, col=coul,las=2)
      }
    }
    else{
      if(j>1 && j<m){
        if(is.null(nr_sub)==F){
          par(mar = c(6, 1, 4, 2))
          barplot(t(a[[j]]),xlab=paste("state",j),yaxt="n",col=coul,las=2)
        }
        else{
          par(mar = c(6, 1, 4, 2))
          barplot(t(a[[j]]),xlab=paste("state",j),yaxt="n", col=coul,las=2)
        }
      }
      else{
        if(is.null(nr_sub)==F){
          par(mar = c(6, 1, 4, 2))
          barplot(t(a[[j]]),xlab=paste("state",j),yaxt="n",legend.text = colnames(a[[j]]),args.legend = list(x="right",inset=-0.2,cex=0.7,bty="n"),col=coul,las=2)
        }
        else{
          par(mar = c(6, 1, 4, 2))
          barplot(t(a[[j]]),xlab=paste("state",j),yaxt="n",legend.text = colnames(a[[j]]),args.legend = list(x="right",inset=-0.2,cex=0.7,bty="n"), col=coul,las=2)
        }
      }

    }


  }
}

barplot_emiss(emiss=subject_emiss,q_emiss = 3,m=3,n_dep =4,group_labs = group,color = coul,nr_sub = 3)
barplot_emiss(emiss=group1,q_emiss = q_emiss1,m=m1,n_dep=n_dep1,group_labs = group,color = coul)

