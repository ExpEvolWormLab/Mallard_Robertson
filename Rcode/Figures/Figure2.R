rm(list=ls());gc()
library(MCMCglmm)

load("output_files/RData/G_matrices_high_salt.RData")
load("output_files/RData/G_matrices_low_salt.RData")
vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")

null_matrix=new.env()
load("output_files/RData_Random/Random_G_Analysis_Cemee_Pop_WI_A6140_NaCl_LIGHT.RData",envir=null_matrix)
A6140_NaCl_null = null_matrix$df_G1
load("output_files/RData_Random/Random_G_Analysis_Cemee_Pop_WI_A6140_NGM_LIGHT.RData",envir=null_matrix)
A6140_NGM_null = null_matrix$df_G1
rm(null_matrix);gc()


######################################################################
##### Eigen decomposition of the NGM/NaCl matrices     ####
######################################################################

sampled_variance_NGM=NULL
sampled_variance_NaCl=NULL

sampled_variance_NGM_null=NULL
sampled_variance_NaCl_null=NULL

sampled_variance_NGM_null_v2=NULL
sampled_variance_NaCl_null_v2=NULL
Lambda_NGM <- eigen(VCV_mat_NGM[[1]]$G1_mat/2)$vectors
Lambda_NaCl <- eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$vectors


for(i in 1:nrow(VCV_mat_NGM[[1]]$VCV_Mat)){
  sampled_variance_NGM =rbind(sampled_variance_NGM,c(sum(eigen(matrix(VCV_mat_NGM[[1]]$VCV_Mat[i,1:49],7,7)/2)$values),eigen(matrix(VCV_mat_NGM[[1]]$VCV_Mat[i,1:49],7,7)/2)$values))
  sampled_variance_NaCl=rbind(sampled_variance_NaCl,c(sum(eigen(matrix(VCV_mat_NaCl[[1]]$VCV_Mat[i,1:49],7,7)/2)$values),eigen(matrix(VCV_mat_NaCl[[1]]$VCV_Mat[i,1:49],7,7)/2)$values))
  
  if(i<=nrow(A6140_NGM_null)){
    sampled_variance_NGM_null =rbind(sampled_variance_NGM_null,c(sum(eigen(matrix(A6140_NGM_null[i,],7,7)/2)$values),eigen(matrix(A6140_NGM_null[i,],7,7)/2)$values))
    sampled_variance_NaCl_null=rbind(sampled_variance_NaCl_null,c(sum(eigen(matrix(A6140_NaCl_null[i,],7,7)/2)$values),eigen(matrix(A6140_NaCl_null[i,],7,7)/2)$values))
    
    ## Rotate the matrix along the true eigenvectors
    sampled_variance_NGM_null_v2  = rbind(sampled_variance_NGM_null_v2,  diag(t(Lambda_NGM)  %*% (matrix(A6140_NGM_null[i,],7,7)/2)  %*% Lambda_NGM))
    sampled_variance_NaCl_null_v2 = rbind(sampled_variance_NaCl_null_v2, diag(t(Lambda_NaCl) %*% (matrix(A6140_NaCl_null[i,],7,7)/2) %*% Lambda_NaCl))
    
  }
}

#### all code have been run to produce the Figures

plot_G_coefs <- function(true_G_NGM, true_G_NaCl,pop_name){
  
  par(mar=c(5,7,4,2))
  vect_Var <- c(2:7,10:14,18:21,26:28,34,35,42,1,9,17,25,33,41,49)
  vProb <- .95
  
  vectX_NGM <- c(true_G_NGM$G1_mat/2)[vect_Var]
  vectX_NaCl <- c(true_G_NaCl$G1_mat/2)[vect_Var]
  
  plot(vectX_NGM,c(31:11,7:1),yaxt="n",bty="n",xlim=c(-.2,.25),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)
  
  lines(c(0,0),c(32.5,9.5))
  lines(c(0,0),c(.5,7.5),col="black")
  
  axis(side=1,pos=0)
  
  axis(side=2,at=c(31:11,7:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF","SF*Size",
                                       "SB*FS","SB*FB","SB*BS","SB*BF","SB*Size",
                                       "FS*FB","FS*BS","FS*BF","FS*Size",
                                       "FB*BS","FB*BF","FB*Size","BS*BF","BS*Size","BF*Size",
                                       "SF","SB","FS","FB","BS","BF","Size"),las=1)
  
  temp_95 <- HPDinterval(true_G_NGM$VCV_Mat[,1:49]/2,prob=.95)
  arrows(temp_95[vect_Var,1],c(31:11,7:1),temp_95[vect_Var,2],c(31:11,7:1),code=3,length=.03,angle=90,lwd=2)
  
  temp_95 <- HPDinterval(true_G_NaCl$VCV_Mat[,1:49]/2,prob=.95)
  arrows(temp_95[vect_Var,1],c(31:11,7:1)+.3,temp_95[vect_Var,2],c(31:11,7:1)+.3,code=3,length=.03,angle=90,lwd=2)
  
  temp_80 <- HPDinterval(true_G_NGM$VCV_Mat[,1:49]/2,prob=.83)
  arrows(temp_80[vect_Var,1],c(31:11,7:1), temp_80[vect_Var,2],c(31:11,7:1),code=3,length=0,angle=90,lwd=3,col="firebrick3")
  
  temp_80 <- HPDinterval(true_G_NaCl$VCV_Mat[,1:49]/2,prob=.83)
  arrows(temp_80[vect_Var,1],c(31:11,7:1)+.3, temp_80[vect_Var,2],c(31:11,7:1)+.3,code=3,length=0,angle=90,lwd=3,col="grey")
  
  points(vectX_NGM,c(31:11,7:1),pch=21,bg="black",cex=.6)
  points(vectX_NaCl,c(31:11,7:1)+.3,pch=21,bg="black",cex=.6)
  
  
  legend(.06,12,c("Low Salt","High Salt"),lwd=2,col=c("firebrick3","grey"))#,pt.bg="grey",pch=c(8,21))
  
}



##### Now the Figure 2

pdf(file='plots/Figure2.pdf',h=8,w=8)

layout(mat=matrix(c(1,1,2,3),2,2), widths = c(1.8,1.5),heights=c(5,6))

par(mar=c(5,4,4,0))
plot_G_coefs(VCV_mat_NGM[[1]],VCV_mat_NaCl[[1]],"A6140")


### Second and third plot
plot(c(0,0.2),c(
  sum(eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values),
  sum(eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values)),bg=c("black","black"),pch=21,ylab='Total genetic variance',xlim=c(-.2,.4),ylim=c(0,.83),type="n",bty="n",xaxt="n",xlab="")

axis(1,at=c(0,.2),c("Low salt","High salt"),las=1)
temp_int <- HPDinterval(mcmc(sampled_variance_NGM[,1]),prob=.95)
arrows(0,temp_int[1],0,temp_int[2],length=.05,code=3,angle=90)
temp_int <- HPDinterval(mcmc(sampled_variance_NaCl[,1]),prob=.95)
arrows(0.2,temp_int[1],0.2,temp_int[2],length=.05,code=3,angle=90)

temp_int <- HPDinterval(mcmc(sampled_variance_NGM[,1]),prob=.83)
arrows(0,temp_int[1],0,temp_int[2],length=0,code=3,angle=90,lwd=2,col="firebrick3")
temp_int <- HPDinterval(mcmc(sampled_variance_NaCl[,1]),prob=.83)
arrows(0.2,temp_int[1],0.2,temp_int[2],length=0,code=3,angle=90,lwd=2,col="grey")

points(c(0,0.2),c(
  sum(eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values),
  sum(eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values)),bg=c("black","black"),pch=21)

###########

plot(eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values~c(1:7),pch=16,xlim=c(1,7.5),bty="n",ylim=c(0,.83),ylab=expression(paste('Genetic variance (',lambda[i],")")),type="n",xaxt="n",xlab="")
axis(1,at=1:7,c(expression(g[max]),expression(g[2]),expression(g[3]),expression(g[4]),expression(g[5]),expression(g[6]),expression(g[7])),las=1)

for(i in 1:7){
  temp_int <- HPDinterval(mcmc(sampled_variance_NGM[,(i+1)]),prob=.95)
  arrows(i,temp_int[1],i,temp_int[2],length=.05,code=3,angle=90)
  temp_int <- HPDinterval(mcmc(sampled_variance_NaCl[,(i+1)]),prob=.95)
  arrows(i+0.2,temp_int[1],i+0.2,temp_int[2],length=.05,code=3,angle=90)
  
  temp_int <- HPDinterval(mcmc(sampled_variance_NGM[,(i+1)]),prob=.83)
  arrows(i,temp_int[1],i,temp_int[2],length=0,code=3,angle=90,lwd=2,col="firebrick3")
  temp_int <- HPDinterval(mcmc(sampled_variance_NaCl[,(i+1)]),prob=.83)
  arrows(i+0.2,temp_int[1],i+0.2,temp_int[2],length=0,code=3,angle=90,lwd=2,col="grey")
  
}

points( (c(1:7)),eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values ,pch=16,col="black")
points( (c(1:7)+0.2),eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values ,pch=16,col="black")

dev.off()

pdf("plots/Figure2_figure_supplement2.pdf",w=5)

temp_95_NaCl <- HPDinterval(VCV_mat_NaCl[[1]]$VCV_Mat[,1:49]/2,prob=.95)
temp_95_NGM <- HPDinterval(VCV_mat_NGM[[1]]$VCV_Mat[,1:49]/2,prob=.95)
temp_95_NULL_NaCl <- HPDinterval(as.mcmc(A6140_NaCl_null/2),prob=.95)
temp_95_NULL_NGM <- HPDinterval(as.mcmc(A6140_NaCl_null/2),prob=.95)


vect_CoVar <- c(2:7,10:14,18:21,26:28,34,35,42)
vect_Diag <- c(1,9,17,25,33,41,49)

plot(1:7,c(VCV_mat_NaCl[[1]]$G1_mat/2)[vect_Diag],pch=21,col="black",las=1,bty="n",ylab=c("Genetic Variance"),bg='firebrick3',ylim=c(0,.2),type="n",xlim=c(0.7,7),xaxt="n",xlab="Transition rates")
axis(side=1,at=1:7,labels=c("SF","SB","FS","FB","BS","BF","Size"))

arrows(1:7+.15,temp_95_NULL_NaCl[vect_Diag,1],1:7+.15,temp_95_NULL_NaCl[vect_Diag,2],code=3,length=.05,col="orange",angle=90,lwd=2)
points(1:7+.15,c(VCV_mat_NaCl[[1]]$G1_mat/2)[vect_Diag],pch=21,col="black",bg="gray")

arrows(1:7-.15,temp_95_NULL_NGM[vect_Diag,1],1:7-.15,temp_95_NULL_NGM[vect_Diag,2],code=3,length=.05,col="orange",angle=90,lwd=2)
points(1:7-.15,c(VCV_mat_NGM[[1]]$G1_mat/2)[vect_Diag],pch=16,col="firebrick3")

legend(.9,.21,"Posterior mode",bty="n")
legend(1.15,.2,c("Low Salt","High Salt"),bty="n",ncol=2,pch=c(16,21),col=c("firebrick3","black"),pt.bg=c("firebrick3","gray"))
legend(.9,.19,lwd=2,"95% CI of null posterior\nmodes",bty="n",col="orange")
dev.off()


pdf(file='plots/Figure2_figure_supplement3.pdf')

plot(eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values~c(1:7),pch=16,xlim=c(1,7.5),bty="n",ylim=c(0.001,.83),ylab=expression(paste('Genetic variance (',lambda[i],")")),type="n",xaxt="n",xlab="",log="y")
axis(1,at=1:7,c(expression(g[max]),expression(g[2]),expression(g[3]),expression(g[4]),expression(g[5]),expression(g[6]),expression(g[7])),las=1)

temp_int <- HPDinterval(mcmc(sampled_variance_NGM_null[,2:8]),prob=.95)
arrows((1:7),temp_int[,1],(1:7),temp_int[,2],length=.05,code=3,angle=90,col="orange",cex=1.2)
temp_int <- HPDinterval(mcmc(sampled_variance_NaCl_null[,2:8]),prob=.95)
arrows((1:7)+0.2,temp_int[,1],(1:7)+0.2,temp_int[,2],length=.05,code=3,angle=90,col="orange")

temp_int <- HPDinterval(mcmc(sampled_variance_NGM_null_v2),prob=.95)
arrows((1:7),temp_int[,1],(1:7),temp_int[,2],length=.05,code=3,angle=90,col="darkgreen",cex=1.2)
temp_int <- HPDinterval(mcmc(sampled_variance_NaCl_null_v2),prob=.95)
arrows((1:7)+0.2,temp_int[,1],(1:7)+0.2,temp_int[,2],length=.05,code=3,angle=90,col="darkgreen")

points( (c(1:7)),eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values ,pch=16,col="firebrick3")
points( (c(1:7)+0.2),eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values ,pch=21,bg="grey")

legend(3,.9,"Posterior mode",bty="n")
legend(3,.6,c("Low Salt","High Salt"),bty="n",ncol=2,pch=c(16,21),pt.bg=c("firebrick3","grey"),col=c("firebrick3","black"))
legend(3,.3,lwd=2,"95% CI of null posterior\nmode",bty="n",col="orange")


dev.off()

