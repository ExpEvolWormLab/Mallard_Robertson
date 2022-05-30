rm(list=ls());gc()
library(MCMCglmm)
library(psych)
library(ggplot2)
library(dplyr)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(corrplot)
library(Rmisc)
library(nlme)
library(lme4)

## Here we only have data for the lines in NaCl - no population measurements
final_merged =read.table("data/Final_merged_data_withNaCl.txt",h=TRUE,sep="\t")
final_merged_NaCl=subset(final_merged, env_label=="NaCl")

#### We only select the relevant data for this project - no Paris phenotyping
final_merged_NaCl$population=substring(final_merged_NaCl$pop_label,1,5)
table((final_merged_NaCl[,c("population","location_label")]))

final_merged_NaCl$is_2012 = substring(final_merged_NaCl$date_str,1,4)!="2012"	
final_merged_NaCl= subset(final_merged_NaCl, ! (population=='A6140' & location_label=='Paris'))
dim(final_merged_NaCl) # 1039 lines

## Populations have been phenotyped in one single location each - not affected
## when estimating the G-matrices. A correction will be needed when estimating the
## comparison CA50 with GA50 if needed.

### For the G calculation, we simply need to remove the additinal data done in Paris
dim(final_merged_NaCl)
data_for_G = final_merged_NaCl
dim(data_for_G) # 1039 - OK

plot.acfs <- function(x) {
  n <- dim(x)[2]
  par(mfrow = c(ceiling(n/2), 2), mar = c(3, 2, 3, 0.5))
  for (i in 1:n) {
    acf(x[, i], lag.max = 100, main = colnames(x)[i], ci.type = "ma", ylim = c(-0.15, 0.15))
    grid()
  }
}


vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32")

VCV_mat_NaCl = NULL
nb_trait = length(vect_P_traits)
k=0
vect_populations=sort(unique(data_for_G$population))
for (i in vect_populations) {
  print(i)
  temp_final = subset(data_for_G,   population == i)
  
  ### Normalize the env. factors
  for(j in c('temperature',"rel_humidity","logD")){
    temp_final[,j] <- (temp_final[,j]-mean(temp_final[,j]))/sd(temp_final[,j])
  }
  
  
  phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))
  prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
                    R = list(V = phen.var/3, n = nb_trait))
  
  model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~  (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                         rcov = ~us(trait):units, 
                         family = rep("gaussian", nb_trait), data = temp_final, prior = prior_mod, verbose = FALSE,nitt=1500000, burnin=500000,thin=100)
  
  pdf(file = paste0("output_files/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_", 
                    i, ".pdf"), height = 20)
  par(mfrow = c(10, 2), mar = c(2, 2, 1, 0))
  plot(model_MCMC$Sol, auto.layout = F)
  dev.off()
  
  pdf(file = paste0("output_files/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_autocorr_", 
                    i, ".pdf"), height = 10)
  plot.acfs(model_MCMC$Sol)
  dev.off()
  
  post_dist = posterior.mode(model_MCMC$VCV)
  k=k+1
  VCV_mat_NaCl[[k]]=list(Population = i, N_measurement = nrow(temp_final), G1_mat = matrix(post_dist[1:nb_trait^2], 
                                                                                      nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
                    R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)
}

###Create tables
k=0
for(vect_pop_id in vect_populations){
  k=k+1
  write.table(matrix(paste0(round(1000*VCV_mat_NaCl[[k]]$G1_mat)/1000
                            ,"  [",round(1000*matrix(HPDinterval(VCV_mat_NaCl[[k]]$VCV_Mat)[1:36,1],ncol=6))/1000,";",
                            round(1000*matrix(HPDinterval(VCV_mat_NaCl[[k]]$VCV_Mat)[1:36,2],ncol=6))/1000,"]"),ncol=6),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,file=paste0("output_files/txt/", vect_pop_id,".txt"))
  
  write.table(VCV_mat_NaCl[[k]]$G1_mat, quote=FALSE,sep="\t", row.names=FALSE,col.names=FALSE, file=paste0("output_files/txt/", vect_pop_id,".txt"))
  
}


#######
#### We should also compute the G matrix in NGM - all done in Lisbon 2012-2013
final_merged_NGM=subset(final_merged,env_label=="NGM")
final_merged_NGM <- subset(final_merged_NGM,substring(date_str,1,4)%in%c("2012","2013"))
dim(final_merged_NGM) # 673

final_merged_NGM$population=substring(final_merged_NGM$pop_label,1,5)
final_merged_NGM$is_2012 = substring(final_merged_NGM$date_str,1,4)!="2012"	

data_for_G_NaCl=data_for_G
data_for_G=final_merged_NGM
data_for_G_NGM=data_for_G


VCV_mat_NGM = NULL
nb_trait = length(vect_P_traits)
k=0
vect_populations=sort(unique(data_for_G$population))
for (i in vect_populations) {
  print(i)
  temp_final = subset(data_for_G,   population == i)
  
  ### Normalize the env. factors
  for(j in c('temperature',"rel_humidity","logD")){
    temp_final[,j] <- (temp_final[,j]-mean(temp_final[,j]))/sd(temp_final[,j])
  }
  
  phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))
  prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
                    R = list(V = phen.var/3, n = nb_trait))
  
  model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~ (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                         rcov = ~us(trait):units, 
                         family = rep("gaussian", nb_trait), data = temp_final, prior = prior_mod, verbose = FALSE,nitt=1500000, burnin=500000,thin=100)
  
  pdf(file = paste0("output_files/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_", 
                    i, "_inNGM.pdf"), height = 20)
  par(mfrow = c(10, 2), mar = c(2, 2, 1, 0))
  plot(model_MCMC$Sol, auto.layout = F)
  dev.off()
  
  pdf(file = paste0("output_files/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_autocorr_", 
                    i, "_inNGM.pdf"), height = 10)
  plot.acfs(model_MCMC$Sol)
  dev.off()
  
  post_dist = posterior.mode(model_MCMC$VCV)
  k=k+1
  VCV_mat_NGM[[k]]=list(Population = i, N_measurement = nrow(temp_final), G1_mat = matrix(post_dist[1:nb_trait^2], 
                                                                                          nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
                        R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)
}


###Create tables
k=0
for(vect_pop_id in vect_populations){
  k=k+1
  write.table(matrix(paste0(round(1000*VCV_mat_NGM[[k]]$G1_mat)/1000
                            ,"  [",round(1000*matrix(HPDinterval(VCV_mat_NGM[[k]]$VCV_Mat)[1:36,1],ncol=6))/1000,";",
                            round(1000*matrix(HPDinterval(VCV_mat_NGM[[k]]$VCV_Mat)[1:36,2],ncol=6))/1000,"]"),ncol=6),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,file=paste0("output_files/G_mat_tables/", vect_pop_id,"_inNGM.txt"))
  
  write.table(VCV_mat_NGM[[k]]$G1_mat, quote=FALSE,sep="\t", row.names=FALSE,col.names=FALSE, file=paste0("output_files/txt/", vect_pop_id,"_inNGM.txt"))
  
}


save(list=ls(),file="output_files/RData/Analysis_Cemee_Pop_WI_NaCl.RData")






