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


load("output_files/RData/Analysis_Cemee_Pop_WI_NaCl.RData")
load("output_files/RData/Analysis_Cemee_Pop_WI_Line_Means_corrected.RData")

angle_theta <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- 180/pi * as.numeric(acos(dot.prod/(norm.x * norm.y)))
	if(theta>90) theta <- (180-theta)
	as.numeric(theta)
}

gmax_NGM <- eigen(VCV_mat_NGM[[1]]$G1_mat)$vectors[,1]
gmax_NaCl <- eigen(VCV_mat_NaCl[[1]]$G1_mat)$vectors[,1]

delta_P_init <- delta_plasticity_lmer#line_phen_medians[1,7:12]-line_phen_medians[1,1:6]

# We can also generate variation in delta_P vector by sampling with replacement in the lines

final_merged_NaCl$pop_label <- paste(final_merged_NaCl$pop_label,"NaCl",sep="_")
final_merged_NGM$pop_label <- paste(final_merged_NGM$pop_label,"NGM",sep="_")


#selected_col <- c(vect_P_traits,"population", "population2","temperature","rel_humidity","logD","pop_label","date_str","env_label","is_2012")
selected_col <- c(vect_P_traits,"population", "temperature","rel_humidity","logD","pop_label","date_str","env_label","is_2012","location_label")
final_export_all <- rbind(final_merged_NaCl[,selected_col],final_merged_NGM[,selected_col])

final_export_all$pop_label_init <- tstrsplit(final_export_all$pop_label,"_")[[1]]
head(final_export_all)

final_export_all$temperature= (final_export_all$temperature-mean(final_export_all $temperature))/sd(final_export_all $temperature)
final_export_all $rel_humidity= (final_export_all $rel_humidity-mean(final_export_all $rel_humidity))/sd(final_export_all $rel_humidity)
final_export_all $logD= (final_export_all $logD-mean(final_export_all$logD))/sd(final_export_all $logD)

final_export_all_A6140 <- subset(final_export_all,population=="A6140" & location_label=="Lisbon")
dim(final_export_all_A6140)

delta_P_sampled <- NULL
n_iter <- 1000
#idx_rdm_P <- matrix(sample(1:nrow(A6_phen_values_NGM),(n_iter*nrow(A6_phen_values_NGM)),replace=TRUE),n_iter,nrow(A6_phen_values_NGM))

Sampled_delta_plasticity_lmer <- array(data=NA,dim=c(n_iter,6))

for(i in 1:n_iter){
# sample lines with replacement
	sampled_lines <- sample(unique(final_export_all_A6140$pop_label_init),replace=TRUE)
	temp_data <- subset(final_export_all_A6140, pop_label_init%in% sampled_lines)
if(i%%100==0) print(i)
for(j in 1:6){

	temp_mod1 <- lmer(temp_data[, vect_P_traits[j]]~ (temperature+rel_humidity+logD)^3 + env_label + (1|pop_label) +(1|date_str),data= temp_data)
	Sampled_delta_plasticity_lmer[i,j] <- c(- as.numeric(coef(temp_mod1)[[1]][1,][names(coef(temp_mod1)[[1]][1,])=="env_labelNGM"]))

}
}

# delta_P_init differs a little bit with the more robusts sampled delta P means, the other 2 does not differ
delta_P <- colMeans(Sampled_delta_plasticity_lmer)
delta_P_sampled <- Sampled_delta_plasticity_lmer
gmax_proj_NGM <- (t(delta_P)%*%(VCV_mat_NGM[[1]]$G1_mat/2)%*%delta_P)/sum(delta_P^2)
pi_E_NGM <- gmax_proj_NGM/eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values[1] # 18%
pi_0_NGM <- mean(eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values)/eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values[1] # 25%

gmax_proj_NaCl<- (t(delta_P)%*%(VCV_mat_NaCl[[1]]$G1_mat/2)%*%delta_P)/sum(delta_P^2)
pi_E_NaCl <- gmax_proj_NaCl/eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values[1] # 11%
pi_0_NaCl <- mean(eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values)/eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values[1] # 22%

theta_E_NGM  <- angle_theta(delta_P,eigen(VCV_mat_NGM[[1]]$G1_mat/2)$vector[,1]) # 69°
theta_E_NaCl  <- angle_theta(delta_P,eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$vector[,1]) # 76°

vect_rand_piE_NGM <- NULL
vect_rand_piE_NaCl <- NULL
vect_rand_thetaE_NGM <- NULL
vect_rand_thetaE_NaCl <- NULL
vect_rand_pi_0_NGM <- NULL
vect_rand_pi_0_NaCl <- NULL
spld_idx <- sample(1:nrow(VCV_mat_NGM[[1]]$VCV_Mat),n_iter)
k=0
for(i in spld_idx){
k=k+1
	temp_NGM <- matrix(VCV_mat_NGM[[1]]$VCV_Mat[i,1:36],6,6)
	temp_NaCl <- matrix(VCV_mat_NaCl[[1]]$VCV_Mat[i,1:36],6,6)

	temp_gmax_proj <- (t(delta_P_sampled[k,])%*%(temp_NGM/2)%*%delta_P_sampled[k,])/sum(delta_P_sampled[k,]^2)
	vect_rand_piE_NGM <- c(vect_rand_piE_NGM ,temp_gmax_proj/eigen(temp_NGM/2)$values[1])

	temp_gmax_proj <- (t(delta_P_sampled[k,])%*%(temp_NaCl/2)%*%delta_P_sampled[k,])/sum(delta_P_sampled[k,]^2)
	vect_rand_piE_NaCl <- c(vect_rand_piE_NaCl , temp_gmax_proj/eigen(temp_NaCl/2)$values[1])	
	
	vect_rand_thetaE_NGM <- c(vect_rand_thetaE_NGM,angle_theta(delta_P_sampled[k,],eigen(temp_NGM/2)$vector[,1]))
	vect_rand_thetaE_NaCl <- c(vect_rand_thetaE_NaCl,angle_theta(delta_P_sampled[k,],eigen(temp_NaCl/2)$vector[,1]))
	
	vect_rand_pi_0_NGM <- c(vect_rand_pi_0_NGM,mean(eigen(temp_NGM/2)$values)/eigen(temp_NGM/2)$values[1])
	vect_rand_pi_0_NaCl <- c(vect_rand_pi_0_NaCl,mean(eigen(temp_NaCl/2)$values)/eigen(temp_NaCl/2)$values[1])
		
}
rm(temp_NGM, temp_NaCl, temp_gmax_proj)

pdf(file='plots/Pi_Theta_Plasticity_A6140.pdf',width=6, height=4)
layout(t(c(1,2)), width=c(0.6,0.4))

plot(1,1,type="n",xlim=c(0,1.5),ylim=c(0.5,2.5),bty="n",las=1,yaxt="n",xlab=expression(Pi),ylab="",xaxt="n")
axis(1,at=c(0,0.5,1))
mtext(side=3,"A.",at=0,cex=2)
temp_95 <- HPDinterval(as.mcmc(vect_rand_piE_NGM))
temp_80 <- HPDinterval(as.mcmc(vect_rand_piE_NGM),prob=0.8)

arrows(temp_95[1,1],2.3,temp_95[1,2],2.3,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],2.3,temp_80[1,2],2.3,code=3,angle=90,length=0,col="cornflowerblue",lwd=2)
points(mean(vect_rand_piE_NGM),2.3,pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_pi_0_NGM))
temp_80 <- HPDinterval(as.mcmc(vect_rand_pi_0_NGM),prob=0.8)
arrows(temp_95[1,1],2.05,temp_95[1,2],2.05,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],2.05,temp_80[1,2],2.05,code=3,angle=90,length=0,col="cornflowerblue",lwd=2)
points(mean(vect_rand_pi_0_NGM),2.05,pch=8)


temp_95 <- HPDinterval(as.mcmc(vect_rand_piE_NaCl))
temp_80 <- HPDinterval(as.mcmc(vect_rand_piE_NaCl),prob=0.8)

arrows(temp_95[1,1],1,temp_95[1,2],1,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],1,temp_80[1,2],1,code=3,angle=90,length=0,col="firebrick3",lwd=2)
points(mean(vect_rand_piE_NaCl),1,pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_pi_0_NaCl));temp_80 <- HPDinterval(as.mcmc(vect_rand_pi_0_NaCl),prob=0.8)
arrows(temp_95[1,1],.75,temp_95[1,2],.75,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],.75,temp_80[1,2],.75,code=3,angle=90,length=0,col="firebrick3",lwd=2)
points(mean(vect_rand_pi_0_NaCl),.75,pch=8)

mtext(side=2,"Low\n Salt",padj=0,adj=0.75,cex=1.2)
mtext(side=2,"High\n Salt",padj=0,adj=0.15,cex=1.2)

plot(1,1,type="n",xlim=c(0,90),xlab=expression(paste(Theta," (°)")),ylim=c(0.5,2.5),bty="n",las=1,yaxt="n",ylab="",xaxt="n")
axis(1,at=c(0,45,90))
mtext(side=3,"B.",at=0,cex=2)
temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_NGM))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_NGM),prob=0.8)

arrows(temp_95[1,1],2,temp_95[1,2],2,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],2,temp_80[1,2],2,code=3,angle=90,length=0,col="cornflowerblue",lwd=2)
points(mean(vect_rand_thetaE_NGM),2,pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_NaCl))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_NaCl),prob=0.8)

arrows(temp_95[1,1],1,temp_95[1,2],1,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],1,temp_80[1,2],1,code=3,angle=90,length=0,col="firebrick3",lwd=2)
points(mean(vect_rand_thetaE_NaCl),1,pch=16)

lines(c(45,45),c(.75,2.25),lty=2)
dev.off()

save(list=ls(),file='output_files/RData/Pi_Theta_Plasticity_A6140.RData')









