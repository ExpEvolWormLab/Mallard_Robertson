rm(list=ls())
gc()
library(MCMCglmm)
library(matrixStats)

load("output_files/RData/G_matrices_high_salt.RData")
load('output_files/RData/Gqw_High_Salt.RData')
load('output_files/RData/Gqw_Low_Salt.RData')

angle_theta <- function(x, y) {
  dot.prod <- x %*% y
  norm.x <- norm(x, type = "2")
  norm.y <- norm(y, type = "2")
  theta <- 180/pi * as.numeric(acos(dot.prod/(norm.x * norm.y)))
  as.numeric(theta)
}

## Opposite sign to match the direction of evolution
load('output_files/RData/SSCP_eigenvectors.RData')
dmax = EV_div[,1]
gmax = eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$vectors[,1]

## Load the coefficients from the Manova
coef_manova <- read.table('output_files/txt/Manova.coefficients.txt',sep="\t")
coef_manova_GA_evolution <- coef_manova[2:4,]
coef_manova_GA_evolution_NGM <- coef_manova[2:4,] +
  coef_manova[110:112,]  # Add the interaction

avg_NaCl = colMeans(coef_manova_GA_evolution)
avg_NGM = colMeans(coef_manova_GA_evolution_NGM)

# New sampled mean
all_ratios <- NULL; all_ratios_drift <- NULL
all_ratios_NGM <- NULL; all_ratios_drift_NGM <- NULL

for(i in 1:1000){

  delta_Z_sampled <- matrix(VCV_with_w_NaCl$VCV_Mat[i,1:64]/2,8,8)[1:7,8]
  delta_Z_sampled_NGM <- matrix(VCV_with_w_NGM$VCV_Mat[i,1:64]/2,8,8)[1:7,8]

  all_ratios <- rbind(all_ratios,c(

    as.numeric(coef_manova_GA_evolution[1,])/delta_Z_sampled,
    as.numeric(coef_manova_GA_evolution[2,])/delta_Z_sampled,
    as.numeric(coef_manova_GA_evolution[3,])/delta_Z_sampled,
    avg_NaCl/delta_Z_sampled))

  all_ratios_NGM <- rbind(all_ratios_NGM,c(
    as.numeric(coef_manova_GA_evolution_NGM[1,])/delta_Z_sampled_NGM,
    as.numeric(coef_manova_GA_evolution_NGM[2,])/delta_Z_sampled_NGM,
    as.numeric(coef_manova_GA_evolution_NGM[3,])/delta_Z_sampled_NGM,
    avg_NGM/delta_Z_sampled_NGM))

}

avg_gen_change <- mean(colMedians(all_ratios[,8:28])) # 3.5 generations on average

angle_pred_dmax = NULL
angle_pred_gmax = NULL
angle_pred_GA1 = NULL;angle_pred_GA2 = NULL;angle_pred_GA4 = NULL
angle_pred_GA1_NGM = NULL;angle_pred_GA2_NGM = NULL;angle_pred_GA4_NGM = NULL

angle_avg_NaCl = NULL;angle_avg_NGM = NULL

for(i in 1:nrow(VCV_with_w_NaCl$VCV_Mat)){
  temp_deltaZ = matrix(VCV_with_w_NaCl$VCV_Mat[i,1:64]/2,8,8)[1:7,8]
  temp_deltaZ_NGM = matrix(VCV_with_w_NGM$VCV_Mat[i,1:64]/2,8,8)[1:7,8]
  
  angle_pred_dmax= c(angle_pred_dmax,angle_theta(dmax,temp_deltaZ))
  angle_pred_gmax= c(angle_pred_gmax,angle_theta(gmax,temp_deltaZ))
  angle_pred_GA1= c(angle_pred_GA1,angle_theta(as.numeric(coef_manova_GA_evolution[1,]),temp_deltaZ))
  angle_pred_GA2= c(angle_pred_GA2,angle_theta(as.numeric(coef_manova_GA_evolution[2,]),temp_deltaZ))
  angle_pred_GA4= c(angle_pred_GA4,angle_theta(as.numeric(coef_manova_GA_evolution[3,]),temp_deltaZ))
  
  angle_pred_GA1_NGM= c(angle_pred_GA1_NGM,angle_theta(as.numeric(coef_manova_GA_evolution_NGM[1,]),temp_deltaZ_NGM))
  angle_pred_GA2_NGM= c(angle_pred_GA2_NGM,angle_theta(as.numeric(coef_manova_GA_evolution_NGM[2,]),temp_deltaZ_NGM))
  angle_pred_GA4_NGM= c(angle_pred_GA4_NGM,angle_theta(as.numeric(coef_manova_GA_evolution_NGM[3,]),temp_deltaZ_NGM))
  
  angle_avg_NaCl= c(angle_avg_NaCl,angle_theta(avg_NaCl,temp_deltaZ))
  angle_avg_NGM= c(angle_avg_NGM,angle_theta(avg_NGM,temp_deltaZ))
  
  
}

angle_Rdm=NULL
for(i in 1:4500){
  angle_Rdm=c(angle_Rdm,angle_theta(runif(7,min=(-1),max=1),runif(7,min=(-1),max=1)))
}
all_angles=cbind(angle_pred_GA1,angle_pred_GA2,angle_pred_GA4,angle_pred_GA1_NGM,angle_pred_GA2_NGM,angle_pred_GA4_NGM,angle_avg_NaCl,angle_avg_NGM,angle_Rdm)

## Vect that contains the beta P
vect_betasP_withBL_NaCl <- NULL

vect_betasP_withBL_NaCl_GA1 <- NULL
vect_betasP_withBL_NaCl_GA2 <- NULL
vect_betasP_withBL_NaCl_GA3 <- NULL

vect_betas_withBL_NaCl_no_cov_sampling <- NULL
vect_betas_withBL_NaCl_no_cov_sampling_only5 <- NULL

for(k in 1:nrow(VCV_with_w_NaCl$VCV_Mat)){
  temp_Gzw_NaCl <- matrix(VCV_with_w_NaCl$VCV_Mat[k,1:64]/2,8,8)
  
  vect_betas_withBL_NaCl_no_cov_sampling <- rbind(vect_betas_withBL_NaCl_no_cov_sampling ,t(solve(temp_Gzw_NaCl[1:7,1:7]/2)%*%VCV_with_w_NaCl$G1_mat[1:7,8]/2))	
  vect_betas_withBL_NaCl_no_cov_sampling_only5 <- rbind(vect_betas_withBL_NaCl_no_cov_sampling_only5 ,
                                                        t(solve(temp_Gzw_NaCl[c(1,2,4,6,7),c(1,2,4,6,7)]/2)%*%VCV_with_w_NaCl$G1_mat[c(1,2,4,6,7),8]/2))
  
  vect_betasP_withBL_NaCl <- rbind(vect_betasP_withBL_NaCl ,t(solve(temp_Gzw_NaCl[1:7,1:7]/2)%*%avg_NaCl))	
  vect_betasP_withBL_NaCl_GA1 <- rbind(vect_betasP_withBL_NaCl_GA1 ,t(solve(temp_Gzw_NaCl[1:7,1:7]/2)%*%as.numeric(coef_manova_GA_evolution[1,])))
  vect_betasP_withBL_NaCl_GA2 <- rbind(vect_betasP_withBL_NaCl_GA2 ,t(solve(temp_Gzw_NaCl[1:7,1:7]/2)%*%as.numeric(coef_manova_GA_evolution[2,])))
  vect_betasP_withBL_NaCl_GA3 <- rbind(vect_betasP_withBL_NaCl_GA3 ,t(solve(temp_Gzw_NaCl[1:7,1:7]/2)%*%as.numeric(coef_manova_GA_evolution[3,])))
}


v_col = c("cadetblue1", "cornflowerblue", "slateblue2")

pdf(file='plots/Figure9.pdf',h=8,w=8)
layout(matrix(c(1,1,2,3),2,2),h=c(.8,1))

par(mar=c(5,4,4,2))
plot(solve(VCV_with_w_NaCl$G1_mat[1:7,1:7])%*%(VCV_with_w_NaCl$G1_mat[1:7,8]),ylim=c(-6,4),type="n",ylab="Selection gradients",xlab="",xlim=c(0.5,8),las=1,bty="n",xaxt="n")
axis(side=1,at=1:7,c("SF","SB","FS","FB","BS","BF","Size"))
legend(1,-4.5,c(expression(beta),expression(beta[g])),lwd=2,col=c("blue","gray"),cex=1.4)

#####

temp_95 <- apply(vect_betasP_withBL_NaCl_GA1,2,function(x){HPDinterval(mcmc(x))})/avg_gen_change
temp_80 <- apply(vect_betasP_withBL_NaCl_GA1,2,function(x){HPDinterval(mcmc(x),prob=.8)})/avg_gen_change

arrows(c(1:7)+.2,temp_95[1,],c(1:7)+.2, temp_95[2,],code=3,length=.025,angle=90,col='black')
arrows(c(1:7)+.2,temp_80[1,],c(1:7)+.2, temp_80[2,],code=0,length=.0,angle=90,col=v_col[1],lwd=2)

temp_95 <- apply(vect_betasP_withBL_NaCl_GA2,2,function(x){HPDinterval(mcmc(x))})/avg_gen_change
temp_80 <- apply(vect_betasP_withBL_NaCl_GA2,2,function(x){HPDinterval(mcmc(x),prob=.8)})/avg_gen_change

arrows(c(1:7)+.3,temp_95[1,],c(1:7)+.3, temp_95[2,],code=3,length=.025,angle=90,col='black')
arrows(c(1:7)+.3,temp_80[1,],c(1:7)+.3, temp_80[2,],code=0,length=.0,angle=90,col=v_col[2],lwd=2)

temp_95 <- apply(vect_betasP_withBL_NaCl_GA3,2,function(x){HPDinterval(mcmc(x))})/avg_gen_change
temp_80 <- apply(vect_betasP_withBL_NaCl_GA3,2,function(x){HPDinterval(mcmc(x),prob=.8)})/avg_gen_change

arrows(c(1:7)+.4,temp_95[1,],c(1:7)+.4, temp_95[2,],code=3,length=.025,angle=90,col='black')
arrows(c(1:7)+.4,temp_80[1,],c(1:7)+.4, temp_80[2,],code=0,length=.0,angle=90,col=v_col[3],lwd=2)

# No sampling of the trait covariances with fitness 
temp_95 <- apply(vect_betas_withBL_NaCl_no_cov_sampling,2,function(x){HPDinterval(mcmc(x))})
temp_80 <- apply(vect_betas_withBL_NaCl_no_cov_sampling,2,function(x){HPDinterval(mcmc(x),prob=.8)})

arrows(c(1:7)-.2,temp_95[1,],c(1:7)-.2, temp_95[2,],code=3,length=.05,angle=90,col='black')
arrows(c(1:7)-.2,temp_80[1,],c(1:7)-.2, temp_80[2,],code=0,length=.0,angle=90,col=c(rep("grey",6),"grey"),lwd=2)
temp_Post <- apply(vect_betas_withBL_NaCl_no_cov_sampling,2,function(x){posterior.mode(mcmc(x))})
points(1:7-.2, temp_Post,pch=16,col='black')

abline(h=0,lty=2)

mtext(side=3,"A",at=0,cex=1.2)

par(mar=c(5,4,1,2))
plot(rep(45,4),1:4,xlim=c(0,180),ylim=c(0.5,6.5),bty="n",las=1,pch=21,bg="grey",type="n",
     xlab="Angle (°)",yaxt="n",ylab="",xaxt="n")
axis(1,at=c(0,45,90,135,180))
axis(2,at=c(1,2,3,5.5),c("GA150","GA250","GA350","Random"),las=1)

abline(v=90,lty=2)

temp_arrows <- apply(all_angles,2,function(x){ HPDinterval(as.mcmc(x))})
temp_arrows_80 <- apply(all_angles,2,function(x){ HPDinterval(as.mcmc(x),prob=0.83)})
arrows(temp_arrows[1,c(1:6,9)],c(1:3-.1,1:3+.1,5.5),temp_arrows[2,c(1:6,9)],c(1:3-.1,1:3+.1,5.5),angle=90,length=.05,code=3)
arrows(temp_arrows_80[1,c(1:6,9)],c(1:3-.1,1:3+.1,5.5), temp_arrows_80[2,c(1:6,9)],c(1:3-.1,1:3+.1,5.5),angle=90,length=.05,code=0,lwd=2,col=c(rep(v_col,2),"gray"))
points(apply(all_angles[,1:6],2,function(x){posterior.mode(as.mcmc(x))}),c(1:3-.1,1:3+.1),pch=rep(c(21,8),each=3),bg="grey",cex=1.3)
legend(95,3,c("Low Salt","High Salt"),pch=c(8,21))
mtext(side=3,"B",at=0,cex=1.2,padj=1)

par(mar=c(5,4,1,2))
plot(apply(all_ratios,2,function(x){posterior.mode(as.mcmc(x))}),rep(1:7,4),xlim=c(-2,9),ylim=c(0.5,8),bty="n",las=1,pch=21,bg="grey",type="n",
     xlab="Ratio Observed/Predicted divergence",yaxt="n",ylab="")

axis(2,at=c(1:7)+.3,c("SF","SB","FS","FB","BS","BF","Size"))
lines(c(0,0),c(.8,7.5),lty=2)

temp_arrows <- apply(all_ratios,2,function(x){ HPDinterval(as.mcmc(x))})
temp_arrows_80 <- apply(all_ratios,2,function(x){ HPDinterval(as.mcmc(x),prob=0.8)})
arrows(temp_arrows[1,8:28],rep(1:7)+rep(c(0,.1,.2),each=7),temp_arrows[2,8:28],rep(1:7)+rep(c(0,.1,.2),each=7),angle=90,length=.05,code=3,lwd=rep(c(1,1,1,1.5),each=7))
arrows(temp_arrows_80[1,8:28],rep(1:7)+rep(c(0,.1,.2),each=7), temp_arrows_80[2,8:28],rep(1:7)+rep(c(0,.1,.2),each=7),angle=90,length=.05,code=0,lwd=rep(c(2,2,2,3),each=7),col=rep(c(v_col),each=7))

points(apply(all_ratios[,8:28],2,function(x){median(x)}),rep(1:7,3)+rep(c(0,.1,.2),each=7),pch=21,bg="grey")
mtext(side=3,"C",at=-3,cex=1.2)
dev.off()


pdf(file='plots/Figure9_figure_supplement1.pdf',w=5)
v_col = c("cadetblue1", "cornflowerblue", "slateblue2")
plot(solve(VCV_with_w_NaCl$G1_mat[1:7,1:7])%*%(VCV_with_w_NaCl$G1_mat[1:7,8]),ylim=c(-2,2),type="n",ylab="Selection gradients",xlab="",xlim=c(0.5,8),las=1,bty="n",xaxt="n")
axis(side=1,at=1:7,c("SF","SB","FS","FB","BS","BF","Size"))

temp_95 <- apply(vect_betas_withBL_NaCl_no_cov_sampling,2,function(x){HPDinterval(mcmc(x))})
temp_80 <- apply(vect_betas_withBL_NaCl_no_cov_sampling,2,function(x){HPDinterval(mcmc(x),prob=.8)})
arrows(c(1:7)-.2,temp_95[1,],c(1:7)-.2, temp_95[2,],code=3,length=.05,angle=90,col='black')
arrows(c(1:7)-.2,temp_80[1,],c(1:7)-.2, temp_80[2,],code=0,length=.0,angle=90,col=c(rep("grey",6),"grey"),lwd=2)
temp_Post <- apply(vect_betas_withBL_NaCl_no_cov_sampling,2,function(x){posterior.mode(mcmc(x))})
points(1:7-.2, temp_Post,pch=16,col='black')

temp_95 <- apply(vect_betas_withBL_NaCl_no_cov_sampling_only5,2,function(x){HPDinterval(mcmc(x))})
temp_80 <- apply(vect_betas_withBL_NaCl_no_cov_sampling_only5,2,function(x){HPDinterval(mcmc(x),prob=.8)})
arrows(c(1,2,4,6,7)+.2,temp_95[1,],c(1,2,4,6,7)+.2, temp_95[2,],code=3,length=.05,angle=90,col='black')
arrows(c(1,2,4,6,7)+.2,temp_80[1,],c(1,2,4,6,7)+.2, temp_80[2,],code=0,length=.0,angle=90,col=c(rep("firebrick3",6),"firebrick3"),lwd=2)
temp_Post <- apply(vect_betas_withBL_NaCl_no_cov_sampling_only5,2,function(x){posterior.mode(mcmc(x))})
points(c(1,2,4,6,7)+.2, temp_Post,pch=16,col='firebrick3')
abline(h=0,lty=2)

dev.off()

vect_betas_NaCl <- NULL

for(k in 1:nrow(VCV_with_w_NaCl$VCV_Mat)){
  temp_Gzw_NaCl <- matrix(VCV_with_w_NaCl$VCV_Mat[k,1:64]/2,8,8)
  vect_betas_NaCl <- rbind(vect_betas_NaCl ,t(solve(temp_Gzw_NaCl[1:7,1:7]/2)%*%(temp_Gzw_NaCl[1:7,8]/2)))	

}


pdf(file='plots/Figure9_figure_supplement2.pdf',w=5)
plot(solve(VCV_with_w_NaCl$G1_mat[1:7,1:7])%*%(VCV_with_w_NaCl$G1_mat[1:7,8]),ylim=c(-1.5,1.5),type="n",ylab=expression(paste("Selection gradient  ",beta[g])),xlab="",xlim=c(0.5,7),xaxt="n",las=1,bty="n")
axis(side=1,c("SF","SB","FS","FB","BS","BF","Size"),at=1:7)
temp_95 <- apply(vect_betas_NaCl,2,function(x){HPDinterval(mcmc(x))})
temp_80 <- apply(vect_betas_NaCl,2,function(x){HPDinterval(mcmc(x),prob=.8)})

arrows(c(1:7),temp_95[1,],c(1:7), temp_95[2,],code=3,length=.05,angle=90,col='black')
arrows(c(1:7),temp_80[1,],c(1:7), temp_80[2,],code=0,length=.0,angle=90,col="gray",lwd=2)

temp_Post <- apply(vect_betas_NaCl,2,function(x){posterior.mode(mcmc(x))})
points(1:7, temp_Post,pch=16,col='black')

abline(h=0,lty=2)

dev.off()


## Export the data in a table :

output_df = data.frame(cbind(temp_Post,t(temp_80),t(temp_95)))
rownames(output_df)=NULL
output_df$trait=c(c("SF","SB","FS","FB","BS","BF","Size"))
output_df=output_df[,c(6,1:5)]
names(output_df)[2:6]=c("Post.mode","lower_83","upper_83","lower_95","upper_95")
output_df[,2:6]=round(output_df[,2:6],digits=3)
write.table(output_df,file=("output_files/txt/Selection_gradients_bg.txt"),sep='\t',quote=FALSE,row.names=FALSE)

