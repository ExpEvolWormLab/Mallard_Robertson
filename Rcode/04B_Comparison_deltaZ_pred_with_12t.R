rm(list=ls())
gc()
library(MCMCglmm)
library(lme4)
library(data.table)

load('output_files/RData/Analysis_Cemee_Pop_WI_NaCl.RData')

Beta_12t <- new.env()
load(file='Output_files/RData/Beta_estimates_all_12t.RData',envir= Beta_12t)
VCV_with_w_NGM_12t = Beta_12t$VCV_with_w_NGM_12t
rm(Beta_12t);gc()

load("output_files/RData/Analysis_Cemee_Pop_WI_Line_Means_corrected.RData")

angle_theta <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- 180/pi * as.numeric(acos(dot.prod/(norm.x * norm.y)))
	as.numeric(theta)
}


deltaZ_post_mode_NaCl <- Beta_12t$VCV_with_w_NGM_12t$G1_mat[7:12,13]
deltaZ_post_mode_NGM <- Beta_12t$VCV_with_w_NGM_12t$G1_mat[1:6,13]

#### We need sampled delta for NGM
final_export_all_NGM <- subset(final_export_all,env_label=='NGM')

# Load bootstrapped data
load("output_files/RData/Sampled_delta_NGM.RData")

### Generate random means and then realized phen. change

load('output_files/RData/Pi_Theta_Plasticity_Evolution.RData')

# New sampled mean
all_angles <- NULL; all_angles_drift <- NULL
all_angles_NGM <- NULL; all_angles_drift_NGM <- NULL

delta_Z_sampled_all <- NULL
delta_Z_sampled_NGM_all <- NULL
	for(i in 1:1000){

	delta_Z_sampled <- matrix(VCV_with_w_NGM_12t$VCV_Mat[sample(1:nrow(VCV_with_w_NGM_12t$VCV_Mat),1),1:169],13,13)[7:12,13]
	
	delta_Z_sampled_all=c(delta_Z_sampled_all,delta_Z_sampled)
	all_angles <- rbind(all_angles,c( angle_theta(delta_Z_sampled, Sampled_delta_GA1_lmer[i,]),
	   angle_theta(delta_Z_sampled, Sampled_delta_GA2_lmer[i,]),
	   angle_theta(delta_Z_sampled, Sampled_delta_GA4_lmer[i,]),angle_theta(delta_Z_sampled,Sampled_delta_GA_lmer[i,])))

	delta_Z_sampled_NGM <- matrix(VCV_with_w_NGM_12t$VCV_Mat[sample(1:nrow(VCV_with_w_NGM_12t$VCV_Mat),1),1:169],13,13)[1:6,13]
	delta_Z_sampled_NGM_all=c(delta_Z_sampled_NGM_all,delta_Z_sampled_NGM)
	all_angles_NGM <- rbind(all_angles_NGM,c( angle_theta(delta_Z_sampled_NGM, Sampled_delta_GA1_lmer_NGM[i,]),
	                                  angle_theta(delta_Z_sampled_NGM, Sampled_delta_GA2_lmer_NGM[i,]),
	                                  angle_theta(delta_Z_sampled_NGM, Sampled_delta_GA4_lmer_NGM[i,]),angle_theta(delta_Z_sampled,Sampled_delta_GA_lmer_NGM[i,])))

}

pdf(file='plots/Angles_delta_Z_both_with12t.pdf',h=7,w=4.5)

plot(apply(all_angles,2,function(x){posterior.mode(as.mcmc(x))}),1:4,xlim=c(0,180),ylim=c(0.5,4.5),bty="n",las=1,pch=21,bg="grey",type="n",
     xlab="Angle (°)",yaxt="n",ylab="",xaxt="n")
axis(1,at=c(0,45,90,135,180))

axis(2,at=c(1:4),c("GA150","GA250","GA350","GA - all"),las=1)

abline(v=90,lty=2)

temp_arrows <- apply(all_angles,2,function(x){ HPDinterval(as.mcmc(x))})
temp_arrows_80 <- apply(all_angles,2,function(x){ HPDinterval(as.mcmc(x),prob=0.8)})
arrows(temp_arrows[1,],1:4,temp_arrows[2,],1:4,angle=90,length=.05,code=3)
arrows(temp_arrows_80[1,],1:4, temp_arrows_80[2,],1:4,angle=90,length=.05,code=0,lwd=2,col="firebrick3")
points(apply(all_angles,2,function(x){posterior.mode(as.mcmc(x))}),1:4,pch=21,bg="grey")

temp_arrows <- apply(all_angles_NGM,2,function(x){ HPDinterval(as.mcmc(x))})
temp_arrows_80 <- apply(all_angles_NGM,2,function(x){ HPDinterval(as.mcmc(x),prob=0.8)})
arrows(temp_arrows[1,],1:4+.2,temp_arrows[2,],1:4+.2,angle=90,length=.05,code=3)
arrows(temp_arrows_80[1,],1:4+.2, temp_arrows_80[2,],1:4+.2,angle=90,length=.05,code=0,lwd=2,col="cornflowerblue")
points(apply(all_angles_NGM,2,function(x){posterior.mode(as.mcmc(x))}),1:4+.2,pch=21,bg="grey")

dev.off()

#############################################
###### Same thing with the strength ??  #####
#############################################

deltaZ_post_mode_NaCl <- VCV_with_w_NGM_12t$G1_mat[7:12,13]
deltaZ_post_mode_NGM <- VCV_with_w_NGM_12t$G1_mat[1:6,13]


### Generate random means and then realized phen. change

# New sampled mean
all_ratios <- NULL; all_ratios_drift <- NULL
all_ratios_NGM <- NULL; all_ratios_drift_NGM <- NULL

for(i in 1:1000){

delta_Z_sampled <- matrix(VCV_with_w_NGM_12t$VCV_Mat[sample(1:nrow(VCV_with_w_NGM_12t$VCV_Mat),1),1:169],13,13)[7:12,13]
delta_Z_sampled_NGM <- matrix(VCV_with_w_NGM_12t$VCV_Mat[sample(1:nrow(VCV_with_w_NGM_12t$VCV_Mat),1),1:169],13,13)[1:6,13]

all_ratios <- rbind(all_ratios,c(
	Sampled_delta_GA1_lmer[i,]/delta_Z_sampled,
	Sampled_delta_GA2_lmer[i,]/delta_Z_sampled,
	Sampled_delta_GA4_lmer[i,]/delta_Z_sampled,
	Sampled_delta_GA_lmer[i,]/delta_Z_sampled))

all_ratios_NGM <- rbind(all_ratios_NGM,c(
  Sampled_delta_GA1_lmer_NGM[i,]/delta_Z_sampled_NGM,
  Sampled_delta_GA2_lmer_NGM[i,]/delta_Z_sampled_NGM,
  Sampled_delta_GA4_lmer_NGM[i,]/delta_Z_sampled_NGM,
  Sampled_delta_GA_lmer_NGM[i,]/delta_Z_sampled_NGM))

}

pdf(file="plots/Ratio_observed_predicted_trait_changes_with12t.pdf")

plot(apply(all_ratios,2,function(x){posterior.mode(as.mcmc(x))}),rep(1:6,4),xlim=c(-2,5),ylim=c(0.5,7),bty="n",las=1,pch=21,bg="grey",type="n",
xlab="Ratio Observed/Predicted divergence",yaxt="n",ylab="")

axis(2,at=c(1:6)+.3,c("SF","SB","FS","FB","BS","BF"))
abline(v=0,lty=2)
#for(k in 5:7) points(angle_theta(deltaZ_post_mode_NaCl,line_phen_medians[k,7:12]-line_phen_medians[1,7:12]),k-4,pch=8)

temp_arrows <- apply(all_ratios,2,function(x){ HPDinterval(as.mcmc(x))})
temp_arrows_80 <- apply(all_ratios,2,function(x){ HPDinterval(as.mcmc(x),prob=0.8)})
arrows(temp_arrows[1,],rep(1:6)+rep(c(0,.1,.2,.6),each=6),temp_arrows[2,],rep(1:6)+rep(c(0,.1,.2,.6),each=6),angle=90,length=.05,code=3,lwd=rep(c(1,1,1,1.5),each=6))
arrows(temp_arrows_80[1,],rep(1:6)+rep(c(0,.1,.2,.6),each=6), temp_arrows_80[2,],rep(1:6)+rep(c(0,.1,.2,.6),each=6),angle=90,length=.05,code=0,lwd=rep(c(2,2,2,3),each=6),col=rep(c(v_col[5:7],'red'),each=6))

points(apply(all_ratios,2,function(x){median(x)}),rep(1:6,4)+rep(c(0,.1,.2,.6),each=6),pch=21,bg="grey")
legend(-2,7,c("GA-all","GA150","GA250","GA450"),col=c(v_col[5:7],'red'),lwd=2,cex=1,bty="n")

dev.off()

pdf('plots/Biplots_predictions_with12t.pdf')
plot(posterior.mode(Sampled_delta_GA_lmer),deltaZ_post_mode_NaCl,pch=16,asp=1,bty="n",type="n",xlim=c(-.5,.3),xlab='Observed divergence',ylab='Predicted divergence')
abline(a=0,b=1)

temp_95 = HPDinterval(as.mcmc(Sampled_delta_GA_lmer))
temp_80 = HPDinterval(as.mcmc(Sampled_delta_GA_lmer),prob=.8)
arrows(temp_95[,1],deltaZ_post_mode_NaCl,temp_95[,2],deltaZ_post_mode_NaCl,angle=90,length=.05,code=3)
arrows(temp_80[,1],deltaZ_post_mode_NaCl,temp_80[,2],deltaZ_post_mode_NaCl,angle=90,length=.05,code=0,lwd=2,col='firebrick3')


temp_95 = HPDinterval(VCV_with_w_NGM_12t$VCV_Mat[,c(91,104,117,130,143,156)])
temp_80 = HPDinterval(VCV_with_w_NGM_12t$VCV_Mat[,c(91,104,117,130,143,156)],prob=.8)

arrows(posterior.mode(Sampled_delta_GA_lmer),temp_95[,1],posterior.mode(Sampled_delta_GA_lmer),temp_95[,2],angle=90,length=.05,code=3)
arrows(posterior.mode(Sampled_delta_GA_lmer),temp_80[,1],posterior.mode(Sampled_delta_GA_lmer),temp_80[,2],angle=90,length=.05,code=0,lwd=2,col='firebrick3')

points(posterior.mode(Sampled_delta_GA_lmer),deltaZ_post_mode_NaCl,pch=16)

temp_95 = HPDinterval(as.mcmc(Sampled_delta_GA_lmer_NGM))
temp_80 = HPDinterval(as.mcmc(Sampled_delta_GA_lmer_NGM),prob=.8)
arrows(temp_95[,1],deltaZ_post_mode_NGM,temp_95[,2],deltaZ_post_mode_NGM,angle=90,length=.05,code=3)
arrows(temp_80[,1],deltaZ_post_mode_NGM,temp_80[,2],deltaZ_post_mode_NGM,angle=90,length=.05,code=0,lwd=2,col='cornflowerblue')

temp_95 = HPDinterval(VCV_with_w_NGM_12t$VCV_Mat[,c(13,26,39,52,65,78)])
temp_80 = HPDinterval(VCV_with_w_NGM_12t$VCV_Mat[,c(13,26,39,52,65,78)],prob=.8)

arrows(posterior.mode(Sampled_delta_GA_lmer_NGM),temp_95[,1],posterior.mode(Sampled_delta_GA_lmer_NGM),temp_95[,2],angle=90,length=.05,code=3)
arrows(posterior.mode(Sampled_delta_GA_lmer_NGM),temp_80[,1],posterior.mode(Sampled_delta_GA_lmer_NGM),temp_80[,2],angle=90,length=.05,code=0,lwd=2,col='cornflowerblue')
points(posterior.mode(Sampled_delta_GA_lmer_NGM),deltaZ_post_mode_NGM,pch=16)

abline(h=0,lty=2)

legend(-.45,.3,c("Low Salt","High Salt"),lwd=2,col=c("cornflowerblue","firebrick3"),bty="n")
dev.off()
