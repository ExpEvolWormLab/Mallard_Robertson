rm(list=ls());gc()
library(MCMCglmm)
#library(lme4)
#library(data.table)
#library(matrixStats)


## Additinal code to look at Gzw matrix with body size - not shown in the manuscript

# Load the phenotypic data
load('output_files/RData/Phenotypic_data.RData')

# Load the fertility data
load("data/fertility_with_n.rda")

vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")

mean_relative_fit <- subset(ecoefs, pop=="A6140" & env=="NaCl")
mean_relative_fit$w <- exp(mean_relative_fit$fertility)/mean(exp(mean_relative_fit$fertility))

data_for_G_NGM_A6140 <- subset(final_merged,population=="A6140" & env_label=='NGM')
data_for_G_NaCl_A6140 <- subset(final_merged,population=="A6140" & env_label=='NaCl')

test_unique_NGM = paste(data_for_G_NGM_A6140$pop_label,data_for_G_NGM_A6140$date_str)
test_unique_NaCl = paste(data_for_G_NaCl_A6140$pop_label,data_for_G_NaCl_A6140$date_str)

data_for_G_NGM_A6140  = data_for_G_NGM_A6140[order(test_unique_NGM),]
data_for_G_NaCl_A6140 = data_for_G_NaCl_A6140[order(test_unique_NaCl),]

shared_uniq = test_unique_NGM[test_unique_NGM%in%test_unique_NaCl]

# Separate the data per environment - NGM =low Salt / NaCl = high salt
data_for_G_NGM_A6140 = subset(data_for_G_NGM_A6140, test_unique_NGM%in%test_unique_NaCl)
data_for_G_NaCl_A6140 = subset(data_for_G_NaCl_A6140, test_unique_NaCl%in%test_unique_NGM)

env_covariates = c("temperature","logD","rel_humidity")

### area.F
retained_names=c(vect_P_traits,env_covariates,'date_str',"pop_label","is_2012")

all_data_NaCl <- data_for_G_NaCl_A6140[,names(data_for_G_NaCl_A6140)%in%retained_names]
all_data_NGM <- data_for_G_NGM_A6140[,names(data_for_G_NGM_A6140)%in%retained_names]

mean_relative_fit <- mean_relative_fit[,c("line","w")]
names(mean_relative_fit)[1] = 'pop_label'

# Merge with fecundity
all_data_NaCl_with_w <- (merge(all_data_NaCl, mean_relative_fit))
dim(all_data_NaCl_with_w) # 337 14
all_data_NGM_with_w <- (merge(all_data_NGM, mean_relative_fit))
dim(all_data_NGM_with_w) # 337 14

# Add the fitness to the list of traits for MCMCglmm
vect_P_traits_all=c(vect_P_traits,"w")
nb_trait = length(vect_P_traits_all)


# Compute the Gqw matrices for A6140 - High Salt
k=0

	phen.var = diag(nb_trait) * diag(var(subset(all_data_NaCl_with_w, select = vect_P_traits_all)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC_NaCl <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32,area.F,w)) ~ 
	                          (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = all_data_NaCl_with_w, prior = prior_mod, verbose = TRUE,nitt=1000000, burnin=100000,thin=200)

	post_dist = posterior.mode(model_MCMC_NaCl$VCV)

	VCV_with_w_NaCl <- list(Population = "A6140_NaCl", N_measurement = nrow(all_data_NaCl_with_w), G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC_NaCl$VCV)

	# Save the list - low salt
	save(list=c("VCV_with_w_NaCl"),file='output_files/RData/Gqw_High_Salt.RData')
	
	# Compute the Gqw matrices for A6140 - Low Salt
	
	phen.var = diag(nb_trait) * diag(var(subset(all_data_NGM_with_w, select = vect_P_traits_all)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
	                  R = list(V = phen.var/3, n = nb_trait))
	
	model_MCMC_NGM <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32,area.F,w)) ~ 
	                              (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
	                            rcov = ~us(trait):units, 
	                            family = rep("gaussian", nb_trait), data = all_data_NGM_with_w, prior = prior_mod, verbose = TRUE,nitt=1000000, burnin=100000,thin=200)
	
	post_dist = posterior.mode(model_MCMC_NGM$VCV)
	
	VCV_with_w_NGM <- list(Population = "A6140_NGM", N_measurement = nrow(all_data_NGM_with_w), G1_mat = matrix(post_dist[1:nb_trait^2], 
	                                                                                                                      nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
	                                R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC_NGM$VCV)
	
	# Save the list - low salt
	save(list=c("VCV_with_w_NGM"),file='output_files/RData/Gqw_Low_Salt.RData')
