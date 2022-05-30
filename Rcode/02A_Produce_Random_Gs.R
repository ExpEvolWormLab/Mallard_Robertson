rm(list=ls());gc()
library(MCMCglmm)
library(ggplot2)
library(dplyr)
library(data.table)
library(matrixStats)
library(boot)
library(Rmisc)
library(nlme)
library(parallel)

load("output_files/RData/Analysis_Cemee_Pop_WI_NaCl.RData")

run_parallel_MCMC <- function(list_param){
	i=list_param[[1]]
	temp_final=list_param[[2]]
	nb_trait=6; vect_P_traits <- c("T12", "T13", "T21", "T23", "T31", "T32")
	
	library(MCMCglmm)
	library(dae)
	library(data.table)	

	temp_final$pop_label <- sample(temp_final$pop_label)
	temp_final$date_str <- sample(temp_final$date_str)	

	### Normalize the env. factors
	for(j in c('temperature',"rel_humidity","logD")){
	  temp_final[,j] <- (temp_final[,j]-mean(temp_final[,j]))/sd(temp_final[,j])
	}
	
	phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~ (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = temp_final, prior = prior_mod, verbose = FALSE,nitt=15000, burnin=5000)

	post_dist = posterior.mode(model_MCMC$VCV)

	VCV_mat_temp=list(Population = i, N_measurement = nrow(temp_final), G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)
		return(VCV_mat_temp)
		
}

clust <- makeCluster(25)

##########
# In NaCl #
##########
data_for_G = subset(final_merged_NaCl, ! (population=="A6140" & location_label=='Paris'))
vect_populations=sort(unique(data_for_G$population))

for(k in vect_populations){
  param_list_NaCl=list()
  for(i in 1:1000) param_list_NaCl[[i]] <- list(i=k,temp_final = subset(data_for_G,   population == k))
  
  List_output <-parLapply(clust, param_list_NaCl , run_parallel_MCMC)
  save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_",k,"_NaCl.RData"))
  
  df_G1 <- NULL
  for(i in 1:length(List_output)){
    df_G1 <- rbind(df_G1,c(List_output[[i]]$G1_mat))
  }
  
  rm(List_output);gc()
  save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_",k,"_NaCl_LIGHT.RData"))
  
  
}

rm(param_list_NaCl);gc()
stopCluster(clust)

clust <- makeCluster(25)

##########
# In NGM #
##########
data_for_G=final_merged_NGM
vect_populations=sort(unique(data_for_G$population))

for(k in vect_populations){
param_list_NGM=list()
for(i in 1:1000) param_list_NGM[[i]] <- list(i=k,temp_final = subset(data_for_G,   population == k))

List_output <-parLapply(clust, param_list_NGM , run_parallel_MCMC)
save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_",k,"_NGM.RData"))

df_G1 <- NULL
for(i in 1:length(List_output)){
df_G1 <- rbind(df_G1,c(List_output[[i]]$G1_mat))
}

rm(List_output);gc()
save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_",k,"_NGM_LIGHT.RData"))


}


stopCluster(clust)
