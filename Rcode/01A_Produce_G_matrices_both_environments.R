rm(list=ls());gc()
library(MCMCglmm)
library(matrixStats)


## Here we only have data for the lines in NaCl - no population measurements
final_merged =read.table("data/Final_merged_data_withNaCl.txt",h=TRUE,sep="\t")

### For the G calculation, we simply need to remove the additional data done in Paris
### these were made to allow comparison with another data set and are not relevant here (small n)
dim(final_merged)
data_for_G = final_merged
dim(data_for_G) # 2883

#### We only select the relevant data for this project - only inbred lines (ie. no Paris phenotyping)
final_merged$population=substring(final_merged$pop_label,1,5)

final_merged$is_2012 = substring(final_merged$date_str,1,4)!="2012"	
final_merged= subset(final_merged, !location_label=='Paris' & !location_label=='Paris2')
dim(final_merged) # 1366 measurements

final_merged$population=as.factor(as.character(final_merged$population))

## Populations have been phenotyped in one single location each - not affected
## when estimating the G-matrices. A correction will be needed when estimating the
## comparison CA50 with GA50 if needed. Here CA50 were removed

### Finally, we need to load the body area data
body_area=read.table("data/body_area_means.txt",h=TRUE,sep='\t')
shared_names=names(final_merged)[names(final_merged)%in%names(body_area)]

temp.merged <- (merge(final_merged,body_area,all.x=TRUE))
temp.merged.missing <- subset(temp.merged,is.na(area.F))[,c("exper_data_id",shared_names)]

dim(final_merged)
final_merged=merge(final_merged,body_area)
dim(final_merged) # 1366

final_merged$line_env=paste0(final_merged$pop_label,final_merged$env_label)

# Export tables with sample sizes
write.table(table(final_merged[,c('population','env_label')]),"output_files/txt/Number_measurements.txt",sep='\t',quote=FALSE)
write.table(table(unique(final_merged[,c("population","line_env","env_label")])[,c('population','env_label')]),"output_files/txt/Number_measurements_uniq_lines.txt",sep='\t',quote=FALSE)

vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")

### Finally we scale the body area to match the sd of the transition rates
### currently in mm^2

mean(colSds(as.matrix(final_merged[final_merged$env_label=="NaCl",vect_P_traits[1:6]]))) # 0.71
mean(colSds(as.matrix(final_merged[final_merged$env_label=="NGM",vect_P_traits[1:6]]))) # 0.54
tapply(final_merged$area.F,final_merged$env_label,sd) # 0.011 / 0.016 in NGM - High Salt

mean(c(mean(colSds(as.matrix(final_merged[final_merged$env_label=="NaCl",vect_P_traits[1:6]])))/tapply(final_merged$area.F,final_merged$env_label,sd)[1],
  mean(colSds(as.matrix(final_merged[final_merged$env_label=="NGM",vect_P_traits[1:6]])))/tapply(final_merged$area.F,final_merged$env_label,sd)[2]))

# Mean factor in each environment ~50 

final_merged$area.F = final_merged$area.F * 50

data_for_G_NaCl <- subset(final_merged,env_label=='NaCl')
data_for_G_NGM <- subset(final_merged,env_label=='NGM')

# Output a sample size table for publication

write.table(cbind(c("","","A6140","GA150","GA250","GA450"),rbind(c("Number of phenotyped RILs per population in each environment",""),c("High Salt","Low Salt"),table(unique(final_merged[,c("population","env_label","pop_label")])[,1:2])))
            ,file='output_files/txt/sample_sizes_uniq_lines.txt',quote=FALSE,col.names=FALSE,row.names=FALSE)

# A function to plot the MCMCglmm outputs
plot.acfs <- function(x) {
  n <- dim(x)[2]
  par(mfrow = c(ceiling(n/2), 2), mar = c(3, 2, 3, 0.5))
  for (i in 1:n) {
    acf(x[, i], lag.max = 100, main = colnames(x)[i], ci.type = "ma", ylim = c(-0.15, 0.15))
    grid()
  }
}

## We will now compute the G matrices for each population separately
VCV_mat_NaCl = NULL
nb_trait = length(vect_P_traits)
k=0
vect_populations=c("A6140","GA150","GA250","GA450")
for (i in vect_populations) {
  print(i)
  temp_final = subset(data_for_G_NaCl,   population == i)
  
  ### Normalize the env. factors
  for(j in c('temperature',"rel_humidity","logD")){
    temp_final[,j] <- (temp_final[,j]-mean(temp_final[,j]))/sd(temp_final[,j])
  }
  
  ## Prior
  phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))
  prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
                    R = list(V = phen.var/3, n = nb_trait))
  
  ## Compute the model
  model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32,area.F)) ~  (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                         rcov = ~us(trait):units, 
                         family = rep("gaussian", nb_trait), data = temp_final, prior = prior_mod, verbose = FALSE,nitt=1500000, burnin=500000,thin=100)
  
  ## Produce the auto-correlation plots
  pdf(file = paste0("output_files/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_", 
                    i, ".pdf"), height = 20)
  par(mfrow = c(10, 2), mar = c(2, 2, 1, 0))
  plot(model_MCMC$Sol, auto.layout = F)
  dev.off()
  
  pdf(file = paste0("output_files/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_autocorr_", 
                    i, ".pdf"), height = 10)
  plot.acfs(model_MCMC$Sol)
  dev.off()
  
  ## Save the results in a list
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
  write.table(matrix(paste0(round(1000*VCV_mat_NaCl[[k]]$G1_mat/2)/1000
                            ,"  [",round(1000*matrix(HPDinterval(VCV_mat_NaCl[[k]]$VCV_Mat/2)[1:49,1],ncol=7))/1000,";",
                            round(1000*matrix(HPDinterval(VCV_mat_NaCl[[k]]$VCV_Mat/2)[1:49,2],ncol=7))/1000,"]"),ncol=7),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,file=paste0("output_files/G_mat_tables/", vect_pop_id,"_High_Salt.txt"))
  
  
}

#######
#### We should also compute the G matrix in NGM

VCV_mat_NGM = NULL
nb_trait = length(vect_P_traits)
k=0

for (i in vect_populations) {
  print(i)
  temp_final = subset(data_for_G_NGM,   population == i)
  
  ### Normalize the env. factors
  for(j in c('temperature',"rel_humidity","logD")){
    temp_final[,j] <- (temp_final[,j]-mean(temp_final[,j]))/sd(temp_final[,j])
  }
  
  phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))
  prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
                    R = list(V = phen.var/3, n = nb_trait))
  
  model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32,area.F)) ~ (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                         rcov = ~us(trait):units, 
                         family = rep("gaussian", nb_trait), data = temp_final, prior = prior_mod, verbose = FALSE,nitt=1500000, burnin=500000,thin=100)
  
  pdf(file = paste0("output_files/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_", 
                    i, "_lowSalt.pdf"), height = 20)
  par(mfrow = c(10, 2), mar = c(2, 2, 1, 0))
  plot(model_MCMC$Sol, auto.layout = F)
  dev.off()
  
  pdf(file = paste0("output_files/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_autocorr_", 
                    i, "_lowSalt.pdf"), height = 10)
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
  write.table(matrix(paste0(round(1000*VCV_mat_NGM[[k]]$G1_mat/2)/1000
                            ,"  [",round(1000*matrix(HPDinterval(VCV_mat_NGM[[k]]$VCV_Mat/2)[1:49,1],ncol=7))/1000,";",
                            round(1000*matrix(HPDinterval(VCV_mat_NGM[[k]]$VCV_Mat/2)[1:49,2],ncol=7))/1000,"]"),ncol=7),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,file=paste0("output_files/G_mat_tables/", vect_pop_id,"_Low_Salt.txt"))
  
  
}

## We save the two lists separately
save(list=c("VCV_mat_NGM"),file="output_files/RData/G_matrices_low_salt.RData")
save(list=c("VCV_mat_NaCl"),file="output_files/RData/G_matrices_high_salt.RData")

## We save the phenotypic data
save(list=c("final_merged"),file="output_files/RData/Phenotypic_data.RData")






