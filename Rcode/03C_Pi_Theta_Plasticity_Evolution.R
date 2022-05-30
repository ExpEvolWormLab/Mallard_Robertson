rm(list=ls())
library(data.table)
library(lme4)
library(plyr)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(matrixStats)
library(boot) 
library(MCMCglmm)


#### Comparison of the Plasticity and the Phenotypic change during evolution (delta Z observed)


v_col = c("chartreuse","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")

#### End of the file production
load("output_files/RData/Analysis_Cemee_Pop_WI_NaCl.RData")
load("output_files/RData/Analysis_Cemee_Pop_WI_Line_Means_corrected.RData")

angle_theta <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- 180/pi * as.numeric(acos(dot.prod/(norm.x * norm.y)))
#	if(theta>90) theta <- (180-theta)
	as.numeric(theta)
}

delta_P_init <- delta_plasticity_lmer

delta_Z1_init <- line_phen_medians[2,7:12]-line_phen_medians[2,7:12]
delta_Z2_init <- line_phen_medians[3,7:12]-line_phen_medians[3,7:12]
delta_Z3_init <- line_phen_medians[4,7:12]-line_phen_medians[4,7:12]

# We can also generate variation in delta_P vector by sampling with replacement in the lines

A6_phen_values_NGM   <- subset(mean_phen_values,population=="A6140")[,2:7]
A6_phen_values_NaCl <- subset(mean_phen_values,population=="A6140")[,8:13]

GA1_phen_values_NaCl <- subset(mean_phen_values,population=="GA150")[,8:13]
GA2_phen_values_NaCl <- subset(mean_phen_values,population=="GA250")[,8:13]
GA3_phen_values_NaCl <- subset(mean_phen_values,population=="GA450")[,8:13]

delta_P_sampled <- NULL
delta_P_sampled2 <- NULL
delta_Z1_sampled <- NULL
delta_Z2_sampled <- NULL
delta_Z3_sampled <- NULL

n_iter <- 1000

# Load the 1000 bootstrapped plasticity vectors
temp_env <- new.env()
load('output_files/RData/Pi_Theta_Plasticity_A6140.RData',envir= temp_env)
delta_P_sampled <- temp_env$delta_P_sampled
rm(temp_env)


## Bootstrap the divergence vectors

final_export_all_NaCl <- subset(final_export_all,env_label=='NaCl')

Sampled_delta_GA1_lmer <- array(data=NA,dim=c(n_iter,6))
Sampled_delta_GA2_lmer <- array(data=NA,dim=c(n_iter,6))
Sampled_delta_GA4_lmer <- array(data=NA,dim=c(n_iter,6))
Sampled_delta_GA_lmer <- array(data=NA,dim=c(n_iter,6))

for(i in 1:n_iter){
# sample lines with replacement
	sampled_lines <- sample(unique(final_export_all_NaCl$pop_label_init),replace=TRUE)
	temp_data <- subset(final_export_all_NaCl, pop_label_init%in% sampled_lines)
if(i%%100==0) print(i)
for(j in 1:6){

temp_mod1 <- lmer(temp_data[, vect_P_traits[j]]~ (temperature+rel_humidity+logD)^3 + population + (1|pop_label) +(1|date_str),data= temp_data)

Sampled_delta_GA1_lmer[i,j] <- as.numeric(coef(temp_mod1)[[1]][1,][names(coef(temp_mod1)[[1]][1,])%in%c("populationGA150")])
Sampled_delta_GA2_lmer[i,j] <- as.numeric(coef(temp_mod1)[[1]][1,][names(coef(temp_mod1)[[1]][1,])%in%c("populationGA250")])
Sampled_delta_GA4_lmer[i,j] <- as.numeric(coef(temp_mod1)[[1]][1,][names(coef(temp_mod1)[[1]][1,])%in%c("populationGA450")])

temp_mod1 <- lmer(temp_data[, vect_P_traits[j]]~ (temperature+rel_humidity+logD)^3 + population2 + (1|pop_label) +(1|date_str)+(1|population),data= temp_data)

Sampled_delta_GA_lmer[i,j] <- as.numeric(coef(temp_mod1)[[1]][1,][names(coef(temp_mod1)[[1]][1,])%in%c("population2GA50")])

}
}

delta_P <- colMeans(delta_P_sampled)

delta_Z1 <- colMeans(Sampled_delta_GA1_lmer)
delta_Z2 <- colMeans(Sampled_delta_GA2_lmer)
delta_Z3 <- colMeans(Sampled_delta_GA4_lmer)

delta_Z_all <- colMeans(Sampled_delta_GA_lmer)


vect_rand_thetaZ1_NaCl <- NULL
vect_rand_thetaZ2_NaCl <- NULL
vect_rand_thetaZ3_NaCl <- NULL

vect_rand_thetaZ_all_NaCl <- NULL

vect_rand_theta0_NaCl <- NULL

for(i in 1:n_iter){

	vect_rand_thetaZ1_NaCl <- c(vect_rand_thetaZ1_NaCl, angle_theta(delta_P_sampled[i,], Sampled_delta_GA1_lmer[i,]))
	vect_rand_thetaZ2_NaCl <- c(vect_rand_thetaZ2_NaCl, angle_theta(delta_P_sampled[i,], Sampled_delta_GA2_lmer[i,]))
	vect_rand_thetaZ3_NaCl <- c(vect_rand_thetaZ3_NaCl, angle_theta(delta_P_sampled[i,], Sampled_delta_GA4_lmer[i,]))		
	vect_rand_thetaZ_all_NaCl <- c(vect_rand_thetaZ_all_NaCl, angle_theta(delta_P_sampled[i,], Sampled_delta_GA_lmer[i,]))		

	vect_rand_theta0_NaCl <- c(vect_rand_theta0_NaCl, angle_theta(delta_P_sampled[sample(1:n_iter,1),], delta_P_sampled[sample(1:n_iter,1),]))		

}

###### Now the plots
pdf(file='plots/Pi_Theta_Plasticity_Evolution_New.pdf',width=6, height=4)
 
par(mar=c(5,7,4,2))
plot(1,1,type="n",xlim=c(0,180),ylim=c(1,2.5),bty="n",las=1,yaxt="n",xlab=expression(paste(Theta," (°)")),ylab="",xaxt="n")#main=expression(Theta [e])
axis(1,at=c(0,45,90,180))
temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaZ1_NaCl))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaZ1_NaCl),prob=0.8)

arrows(temp_95[1,1],1.5,temp_95[1,2],1.5,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],1.5,temp_80[1,2],1.5,code=3,angle=90,length=0,col= v_col[2],lwd=2)
points(mean(vect_rand_thetaZ1_NaCl),1.5,pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaZ2_NaCl))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaZ2_NaCl),prob=0.8)

arrows(temp_95[1,1],1.25,temp_95[1,2],1.25,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],1.25,temp_80[1,2],1.25,code=3,angle=90,length=0,col= v_col[3],lwd=2)
points(mean(vect_rand_thetaZ2_NaCl),1.25,pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaZ3_NaCl))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaZ3_NaCl),prob=0.8)

arrows(temp_95[1,1],1,temp_95[1,2],1,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],1,temp_80[1,2],1,code=3,angle=90,length=0,col= v_col[4],lwd=2)
points(mean(vect_rand_thetaZ3_NaCl),1,pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaZ_all_NaCl))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaZ_all_NaCl),prob=0.8)

arrows(temp_95[1,1],2.25,temp_95[1,2],2.25,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],2.25,temp_80[1,2],2.25,code=3,angle=90,length=0,col="firebrick3",lwd=2)
points(mean(vect_rand_thetaZ_all_NaCl),2.25,pch=16)


abline(v=90,lty=2)

axis(2,at=c(1.5,1.25,1,2.25),labels=c(
  expression(paste(Delta,"P - GA1")),
  expression(paste(Delta,"P - GA2")),
  expression(paste(Delta,"P - GA4")),
  expression(paste(Delta,"P - GA-all"))
    ),las=1)

dev.off()

save(list=ls(),file='output_files/RData/Pi_Theta_Plasticity_Evolution.RData')






