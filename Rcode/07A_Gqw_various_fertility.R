rm(list=ls())
gc()
#load libraries
library(MCMCglmm)
library(parallel)
# Load the phenotypic data
load('output_files/RData/Phenotypic_data.RData')

# Load the fertility data
load("data/fertility_with_n.rda")

vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")

mean_relative_fit <- subset(ecoefs, pop%in%c("A6140","GA") & env=="NaCl")
names(mean_relative_fit)[4] <- "Std_Error"
mean_relative_fit$Sd <- mean_relative_fit[,"Std_Error"] * sqrt(2) 

# For each line, we need to produce an entire distribution
new_fert_data <- data.frame(pop_label= mean_relative_fit$line)
nline <- length(mean_relative_fit$fertility)
nobs= 1000
new_fert_data = cbind(new_fert_data,matrix(rnorm(n=nobs*nline,mean=mean_relative_fit$fertility,sd=mean_relative_fit$Sd),nline,nobs))
new_fert_data <- new_fert_data[order(new_fert_data$pop_label),]

data_for_G_NGM <- subset(final_merged, env_label=='NGM')
data_for_G_NaCl <- subset(final_merged,env_label=='NaCl')

# We need to center the phenotypic measurements
for(i in vect_P_traits){
  data_for_G_NaCl[,i] <- data_for_G_NaCl[,i] - mean(data_for_G_NaCl[,i])
  data_for_G_NGM[,i] <- data_for_G_NGM[,i] - mean(data_for_G_NGM[,i])
}

# We need to normalize the covariates
for(j in c('temperature',"rel_humidity","logD")){
  data_for_G_NaCl[,j] <- (data_for_G_NaCl[,j]-mean(data_for_G_NaCl[,j]))/sd(data_for_G_NaCl[,j])
  data_for_G_NGM[,j] <- (data_for_G_NGM[,j]-mean(data_for_G_NGM[,j]))/sd(data_for_G_NGM[,j])
}

# Only line for which we have fitness measurement
data_for_Gzw_NaCl= subset(data_for_G_NaCl,population=="A6140" & pop_label%in%mean_relative_fit$line)
dim(data_for_Gzw_NaCl) # 339

temp <- table(data_for_Gzw_NaCl$pop_label)

## How many obervations are needed for each line
nb_obs_per_line <- data.frame(pop_label=names(temp),nb_obs=as.numeric(temp))
nb_obs_per_line = subset(nb_obs_per_line,nb_obs>0)
nb_obs_per_line <- nb_obs_per_line[order(nb_obs_per_line$pop_label),];rm(temp)
sum(nb_obs_per_line$nb_obs) # 339

####
# We will now run 100 estimations of the Gzw matrix with sampled fitness values for each behavioral observation
# and a second one with randomized fitness values accross lines ID
####


#### NaCl first

# And only fitness data for which we have phenotypic measuremnt
new_fert_data_sub = subset(new_fert_data,pop_label%in%data_for_Gzw_NaCl$pop_label)
data_for_Gzw_NaCl=data_for_Gzw_NaCl[order(data_for_Gzw_NaCl$pop_label),]

input_list=list()
for(i in 1:500){
  
  temp_fitness_values <- do.call(rbind, apply(new_fert_data_sub,1,function(x){
    data.frame(pop_label=as.character(x[1]),fertility = sample(as.numeric(x[2:(nobs+1)]),subset(nb_obs_per_line,pop_label==x[1])$nb_obs))
  }))
  
  w = exp(temp_fitness_values$fertility)/mean(exp(temp_fitness_values$fertility))
  
  input_list[[i]] <- data.frame(data_for_Gzw_NaCl,w)
  
}


### Function to run the two models
estimate_Gzw_matrices <- function(input_temp){
  vect_P_traits=c("T12","T13",'T21','T23','T31','T32',"area.F")
  library(MCMCglmm)
  
  nb_trait=length(vect_P_traits)+1
  phen.var = diag(nb_trait) * diag(var(subset(input_temp, select = vect_P_traits)))
  prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
                    R = list(V = phen.var/3, n = nb_trait))
  
  
  model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32,area.F,w)) ~ 
                           (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                         rcov = ~us(trait):units, 
                         family = rep("gaussian", nb_trait), data = input_temp, prior = prior_mod, verbose = TRUE,nitt=150000, burnin=100000,thin=10)
  
  
  input_temp$w <- sample(input_temp$w)
  
    post_dist = posterior.mode(model_MCMC$VCV)
  
    if(length(model_MCMC)==20){
    return(posterior.mode(model_MCMC$VCV))
  }else{
    return(NULL)
  }
  
}

ncores=20
clust <- makeCluster(ncores)

output_list <- parLapply(clust, input_list, estimate_Gzw_matrices)
stopCluster(clust)


save(list=c("output_list"),file="output_files/RData/Gqw_with_various_fertility.RData")

######################## Analyze the results ########################
rm(list=ls())
gc()
library(MCMCglmm)
load("output_files/RData/Gqw_with_various_fertility.RData")
load("output_files/RData/Gqw_High_Salt.RData")
all_Gzw = NULL
for(i in 1:length(output_list)){
  all_Gzw= rbind(all_Gzw,output_list[[i]][1:64])
}

# pdf("plots/Effect_of_fitness_variation.pdf",h=8,w=5)
# par(mar=c(5,7,4,2))
# 
# vect_6t_v2 <- c(c(2:7,11:15,20:23,29:31,38,39,47),c(57:63),c(9*(1:8)-8))
# vProb <- .95
# 
# vectX_VarFit <- c(posterior.mode(as.mcmc(all_Gzw))/2)[vect_6t_v2]
# vectX_NaCl <- c(VCV_with_w_NaCl$G1_mat/2)[vect_6t_v2]
# 
# vect_y_v2 =c(c(11:31),c(1:7),c(35:41),8)
# 
# plot(vectX_VarFit,vect_y_v2,yaxt="n",bty="n",xlim=c(-.25,.25),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2,ylim=c(35,0))
# 
# lines(c(0,0),c(31.5,10.5))
# lines(c(0,0),c(.5,7.5),col="black")
# lines(c(0,0),c(34.5,41.5),col="black")
# 
# 
# axis(side=1,pos=44)
# 
# axis(side=2,at=c(1:8,11:31,35:41),labels=c("SF*w","SB*w","FS*w","FB*w","BS*w","BF*w","Area*w","w",
#                                            "SF*SB","SF*FS","SF*FB","SF*BS","SF*BF","SF*Area",
#                                            "SB*FS","SB*FB","SB*BS","SB*BF","SB*Area",
#                                            "FS*FB","FS*BS","FS*BF","FS*Area",
#                                            "FB*BS","FB*BF","FB*Area","BS*BF","BS*Area","BF*Area",
#                                            "SF","SB","FS","FB","BS","BF","Area"),las=1)
# 
# 
# temp_95 <- HPDinterval(as.mcmc(all_Gzw)/2,prob=.95)
# arrows(temp_95[vect_6t_v2,1],vect_y_v2,temp_95[vect_6t_v2,2],vect_y_v2,code=3,length=.02,angle=90)
# 
# temp_80 <- HPDinterval(as.mcmc(all_Gzw)/2,prob=.83)
# arrows(temp_80[vect_6t_v2,1],vect_y_v2, temp_80[vect_6t_v2,2],vect_y_v2,code=3,length=0,angle=90,lwd=2,col=
#          rep(c("cornflowerblue","darkgreen","cornflowerblue"),c(6,1,21)))
# 
# points(vectX_VarFit,vect_y_v2,pch=21,bg="black",cex=.6)
# points(vectX_NaCl,vect_y_v2+.3,pch=21,bg="black",cex=.6)
# 
# ## Kept from the old one
# legend(-.28,25,c("Various fitness","Mean fitness"),lwd=2,lty=c(1,0),col=c("cornflowerblue","black"),cex=.8,pch=16)
# 
# dev.off()


pdf("plots/Figure4_supplement_figure2.pdf",w=6)
par(mar=c(5,7,4,2))

vect_6t_v2_VarFit <- c(57:64)
vect_6t_v2_NaCl <- c(57:64)
vectX_12t_v2 <- c(posterior.mode(as.mcmc(all_Gzw))[vect_6t_v2_VarFit]/2,VCV_with_w_NaCl$G1_mat[vect_6t_v2_NaCl]/2)

vProb <- .95
vect_y_v2 =c(1:8-.2,1:8+.2)

plot(vectX_12t_v2,vect_y_v2,ylim=c(8,0),yaxt="n",bty="n",xlim=c(-.2,.4),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)
lines(c(0,0),c(.5,7.5),col="black")
axis(side=1,pos=8.5)
axis(side=2,at=c(1:8),labels=c("SF*w","SB*w","FS*w","FB*w","BS*w","BF*w","Size*w","w"),las=1)

temp_95 <- HPDinterval(as.mcmc(all_Gzw)/2,prob=.95)
arrows(temp_95[vect_6t_v2_VarFit,1],vect_y_v2[1:8],temp_95[vect_6t_v2_VarFit,2],vect_y_v2[1:8],code=3,length=.05,angle=90)
temp_95 <- HPDinterval(VCV_with_w_NaCl$VCV_Mat[,1:64]/2,prob=.95)

temp_80 <- HPDinterval(as.mcmc(all_Gzw)/2,prob=.83)
arrows(temp_80[vect_6t_v2_VarFit,1],vect_y_v2[1:8],temp_80[vect_6t_v2_VarFit,2],vect_y_v2[1:8],code=3,length=0,angle=90,lwd=2,col="darkgreen")

temp_80 <- HPDinterval(VCV_with_w_NaCl$VCV_Mat[,1:64]/2,prob=.83)

points(vectX_12t_v2,vect_y_v2,bg="black",pch=rep(c(8,16),each=8))
legend(.1,5,c("Variable self-fertility","Average self-fertility"),lwd=2,lty=c(1,0),col=c("darkgreen","black"),cex=1,pch=c(8,16))
dev.off()





