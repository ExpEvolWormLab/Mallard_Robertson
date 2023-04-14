rm(list=ls());gc()

col_GA <- c("cadetblue1", "cornflowerblue", "slateblue2")
load("output_files/RData/G_matrices_high_salt.RData")

sampled_variance_NaCl=NULL
sampled_variance_GA1=NULL
sampled_variance_GA2=NULL
sampled_variance_GA4=NULL

for(i in 1:nrow(VCV_mat_NaCl[[1]]$VCV_Mat)){
  
  sampled_variance_NaCl=rbind(sampled_variance_NaCl,c(sum(eigen(matrix(VCV_mat_NaCl[[1]]$VCV_Mat[i,1:49],7,7)/2)$values),eigen(matrix(VCV_mat_NaCl[[1]]$VCV_Mat[i,1:49],7,7)/2)$values))
  
  sampled_variance_GA1 =rbind(sampled_variance_GA1,c(sum(eigen(matrix(VCV_mat_NaCl[[2]]$VCV_Mat[i,1:49],7,7)/2)$values),eigen(matrix(VCV_mat_NaCl[[2]]$VCV_Mat[i,1:49],7,7)/2)$values))
  sampled_variance_GA2 =rbind(sampled_variance_GA2,c(sum(eigen(matrix(VCV_mat_NaCl[[3]]$VCV_Mat[i,1:49],7,7)/2)$values),eigen(matrix(VCV_mat_NaCl[[3]]$VCV_Mat[i,1:49],7,7)/2)$values))
  sampled_variance_GA4 =rbind(sampled_variance_GA4,c(sum(eigen(matrix(VCV_mat_NaCl[[4]]$VCV_Mat[i,1:49],7,7)/2)$values),eigen(matrix(VCV_mat_NaCl[[4]]$VCV_Mat[i,1:49],7,7)/2)$values))
  

}

all_sampled_variance=list(sampled_variance_NaCl,sampled_variance_GA1,sampled_variance_GA2,sampled_variance_GA4)

load('output_files/RData/Random_skewers.RData')
pdf(file='plots/Figure7.pdf',h=8,w=7.5)

layout(matrix(c(1,1,2,3),2,2),w=c(1.2,1))

par(mar=c(5,7,4,2))
vect_Var <- c(2:7,10:14,18:21,26:28,34,35,42,1,9,17,25,33,41,49)
vProb <- .95

plot(c(VCV_mat_NaCl[[1]]$G1_mat/2)[vect_Var],c(31:11,7:1),yaxt="n",bty="n",xlim=c(-.25,.50),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)
mtext(side=3,at=-.5,"A",cex=1.5)
lines(c(0,0),c(31.5,8.5))
lines(c(0,0),c(.5,7.5))

axis(side=1,pos=0)

axis(side=2,at=c(31:11,7:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF","SF*Size",
                                     "SB*FS","SB*FB","SB*BS","SB*BF","SB*Size",
                                     "FS*FB","FS*BS","FS*BF","FS*Size",
                                     "FB*BS","FB*BF","FB*Size","BS*BF","BS*Size","BF*Size",
                                     "SF","SB","FS","FB","BS","BF","Size"),las=1)

for(i in 1:4){
  
  temp_95 <- HPDinterval(VCV_mat_NaCl[[i]]$VCV_Mat[,1:49]/2,prob=.95)
  arrows(temp_95[vect_Var,1],c(31:11,7:1)+(.15*(i-1)),temp_95[vect_Var,2],c(31:11,7:1)+(.15*(i-1)),code=3,length=.02,angle=90)
  
  temp_80 <- HPDinterval(VCV_mat_NaCl[[i]]$VCV_Mat[,1:49]/2,prob=.83)
  arrows(temp_80[vect_Var,1],c(31:11,7:1)+(.15*(i-1)),
         temp_80[vect_Var,2],c(31:11,7:1)+(.15*(i-1)),code=3,length=0,angle=90,lwd=2,col=c("grey",col_GA)[i])
  
  points(c(VCV_mat_NaCl[[i]]$G1_mat/2)[vect_Var],c(31:11,7:1)+(.15*(i-1)),pch=21,bg="black",cex=.6)
}

legend(.15,23,c("A6140","GA150","GA250","GA450"),lwd=2,bty="s",col=c("grey",col_GA))

par(mar=c(5,4,4,2))

plot(c(0,1),c(
  sum(eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values),
  sum(eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values)),bg=c("black","black"),pch=21,ylab='Total genetic variance',xlim=c(0,1.3),ylim=c(0,.83),type="n",bty="n",xaxt="n",xlab="")
#axis(side=1,at=c(0,1),labels=c("A6140","GA[1,2,4]50"),padj=.5)
mtext(side=3,at=-.5,"B",cex=1.5)
posX=c(0,.7,1,1.3)
axis(side=1,at=posX,c("A6140","GA150","GA250","GA450"),las=2)

for(i in 1:4){
  temp_int <- HPDinterval(mcmc(all_sampled_variance[[i]][,1]),prob=.95)
  arrows(posX[i],temp_int[1],posX[i],temp_int[2],length=.05,code=3,angle=90)
  temp_int <- HPDinterval(mcmc(all_sampled_variance[[i]][,1]),prob=.83)
  arrows(posX[i],temp_int[1],posX[i],temp_int[2],length=0,code=3,angle=90,lwd=2,col=c("grey",col_GA)[i])
  
}


points(posX,c(
  sum(eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values),
  sum(eigen(VCV_mat_NaCl[[2]]$G1_mat/2)$values),sum(eigen(VCV_mat_NaCl[[3]]$G1_mat/2)$values),sum(eigen(VCV_mat_NaCl[[4]]$G1_mat/2)$values)),bg=c("grey",col_GA),pch=21)

# Variance along e_max

plot(posX,c(0,0,0,2),type="n",las=1,bty="n",ylab=expression(paste("Genetic variance along ",e[max])),xaxt="n",xlab="")
mtext(side=3,at=-.5,"C",cex=1.5)
axis(side=1,at=posX,c("A6140","GA150","GA250","GA450"),las=2)
arrows(posX,HPD.R.vec.proj[,1,1],posX,HPD.R.vec.proj[,2,1],code=3,angle=90,length=0.05)
arrows(posX,HPD.R.vec.proj.83[,1,1],posX,HPD.R.vec.proj.83[,2,1],length=0,col=c("grey",col_GA),lwd=2)
points(posX,R.postmode[1,],pch=16)

dev.off()

# Load the low salt G matrices
load("output_files/RData/G_matrices_low_salt.RData")

## Load the randomized G matrices
null_matrix=new.env()
load("output_files/RData_Random/Random_G_Analysis_Cemee_Pop_WI_A6140_NaCl_LIGHT.RData",envir=null_matrix)
A6140_NaCl_null = null_matrix$df_G1
load("output_files/RData_Random/Random_G_Analysis_Cemee_Pop_WI_A6140_NGM_LIGHT.RData",envir=null_matrix)
A6140_NGM_null = null_matrix$df_G1
rm(null_matrix);gc()

sampled_variance_NGM_null=NULL
sampled_variance_NaCl_null=NULL

sampled_variance_NGM_null_v2=NULL
sampled_variance_NaCl_null_v2=NULL

Lambda_NGM <- eigen(VCV_mat_NGM[[1]]$G1_mat/2)$vectors
Lambda_NaCl <- eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$vectors

for(i in 1:nrow(A6140_NGM_null)){
  sampled_variance_NGM_null =rbind(sampled_variance_NGM_null,c(sum(eigen(matrix(A6140_NGM_null[i,],7,7)/2)$values),eigen(matrix(A6140_NGM_null[i,],7,7)/2)$values))
  sampled_variance_NaCl_null=rbind(sampled_variance_NaCl_null,c(sum(eigen(matrix(A6140_NaCl_null[i,],7,7)/2)$values),eigen(matrix(A6140_NaCl_null[i,],7,7)/2)$values))

  ## Rotate the matrix along the true eigenvectors
  sampled_variance_NGM_null_v2  = rbind(sampled_variance_NGM_null_v2,  diag(t(Lambda_NGM)  %*% (matrix(A6140_NGM_null[i,],7,7)/2)  %*% Lambda_NGM))
  sampled_variance_NaCl_null_v2 = rbind(sampled_variance_NaCl_null_v2, diag(t(Lambda_NaCl) %*% (matrix(A6140_NaCl_null[i,],7,7)/2) %*% Lambda_NaCl))

}


### Figure Supplement
### We need the null G matrices of the GA[1-4] populations

null_matrix=new.env()
load("output_files/RData_Random/Random_G_Analysis_Cemee_Pop_WI_GA150_NaCl_LIGHT.RData",envir=null_matrix)
GA150_NaCl_null = null_matrix$df_G1
load("output_files/RData_Random/Random_G_Analysis_Cemee_Pop_WI_GA250_NaCl_LIGHT.RData",envir=null_matrix)
GA250_NaCl_null = null_matrix$df_G1
load("output_files/RData_Random/Random_G_Analysis_Cemee_Pop_WI_GA450_NaCl_LIGHT.RData",envir=null_matrix)
GA450_NaCl_null = null_matrix$df_G1
rm(null_matrix);gc()

sampled_variance_GA1_null=NULL
sampled_variance_GA2_null=NULL
sampled_variance_GA4_null=NULL

sampled_variance_GA1_null_v2=NULL
sampled_variance_GA2_null_v2=NULL
sampled_variance_GA4_null_v2=NULL

Lambda_GA1 <- eigen(VCV_mat_NaCl[[2]]$G1_mat/2)$vectors
Lambda_GA2 <- eigen(VCV_mat_NaCl[[3]]$G1_mat/2)$vectors
Lambda_GA4 <- eigen(VCV_mat_NaCl[[4]]$G1_mat/2)$vectors

for(i in 1:nrow(GA150_NaCl_null)){
  sampled_variance_GA1_null =rbind(sampled_variance_GA1_null,c(sum(eigen(matrix(GA150_NaCl_null[i,],7,7)/2)$values),eigen(matrix(GA150_NaCl_null[i,],7,7)/2)$values))
  sampled_variance_GA2_null =rbind(sampled_variance_GA2_null,c(sum(eigen(matrix(GA250_NaCl_null[i,],7,7)/2)$values),eigen(matrix(GA250_NaCl_null[i,],7,7)/2)$values))
  sampled_variance_GA4_null =rbind(sampled_variance_GA4_null,c(sum(eigen(matrix(GA450_NaCl_null[i,],7,7)/2)$values),eigen(matrix(GA450_NaCl_null[i,],7,7)/2)$values))
  
  ## Rotate the matrix along the true eigenvectors of A6140
  sampled_variance_GA1_null_v2  = rbind(sampled_variance_GA1_null_v2,  diag(t(Lambda_NaCl)  %*% (matrix(GA150_NaCl_null[i,],7,7)/2)  %*% Lambda_NaCl))
  sampled_variance_GA2_null_v2  = rbind(sampled_variance_GA2_null_v2,  diag(t(Lambda_NaCl)  %*% (matrix(GA250_NaCl_null[i,],7,7)/2)  %*% Lambda_NaCl))  
  sampled_variance_GA4_null_v2  = rbind(sampled_variance_GA4_null_v2,  diag(t(Lambda_NaCl)  %*% (matrix(GA450_NaCl_null[i,],7,7)/2)  %*% Lambda_NaCl))  
}

true_GA1_v2 <- diag(t(Lambda_NaCl)  %*% (VCV_mat_NaCl[[2]]$G1_mat/2)  %*% Lambda_NaCl)
true_GA2_v2 <- diag(t(Lambda_NaCl)  %*% (VCV_mat_NaCl[[3]]$G1_mat/2)  %*% Lambda_NaCl)
true_GA4_v2 <- diag(t(Lambda_NaCl)  %*% (VCV_mat_NaCl[[4]]$G1_mat/2)  %*% Lambda_NaCl)

pdf(file='plots/Figure7_figure_supplement1.pdf')#,h=9,w=4)
temp_95_NULL_GA1 <- HPDinterval(as.mcmc(GA150_NaCl_null/2),prob=.95)
temp_95_NULL_GA2 <- HPDinterval(as.mcmc(GA250_NaCl_null/2),prob=.95)
temp_95_NULL_GA4 <- HPDinterval(as.mcmc(GA450_NaCl_null/2),prob=.95)


vect_CoVar <- c(2:7,10:14,18:21,26:28,34,35,42)
vect_Diag <- c(1,9,17,25,33,41,49)

plot(1:7,c(VCV_mat_NaCl[[2]]$G1_mat/2)[vect_Diag],pch=21,col="black",las=1,bty="n",ylab=c("Genetic Variance"),ylim=c(0,.2),type="n",xlim=c(0.7,7),xaxt="n",xlab="Transition rates")
axis(side=1,at=1:7,labels=c("SF","SB","FS","FB","BS","BF","Size"))

arrows(1:7-.15,temp_95_NULL_GA1[vect_Diag,1],1:7-.15,temp_95_NULL_GA1[vect_Diag,2],code=3,length=.05,col="orange",angle=90,lwd=2)
points(1:7-.15,c(VCV_mat_NaCl[[2]]$G1_mat/2)[vect_Diag],pch=21,col="black",bg=col_GA[1])

arrows(1:7,temp_95_NULL_GA2[vect_Diag,1],1:7,temp_95_NULL_GA2[vect_Diag,2],code=3,length=.05,col="orange",angle=90,lwd=2)
points(1:7,c(VCV_mat_NaCl[[3]]$G1_mat/2)[vect_Diag],pch=21,col="black",bg=col_GA[2])

arrows(1:7+.15,temp_95_NULL_GA4[vect_Diag,1],1:7+.15,temp_95_NULL_GA4[vect_Diag,2],code=3,length=.05,col="orange",angle=90,lwd=2)
points(1:7+.15,c(VCV_mat_NaCl[[4]]$G1_mat/2)[vect_Diag],pch=21,col="black",bg=col_GA[3])


legend(.9,.21,"Posterior mode",bty="n")
legend(1.15,.2,c("GA150","GA250","GA450"),bty="n",ncol=3,pch=21,pt.bg=col_GA,col="black")
legend(.9,.19,lwd=2,"95% CI of null posterior\nmodes",bty="n",col="orange")
dev.off()


pdf(file='plots/Figure7_figure_supplement2.pdf')#,h=9,w=4)
plot(eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values~c(1:7),pch=16,xlim=c(1,7.5),bty="n",ylim=c(0.001,.83),ylab=expression(paste('Genetic variance (',lambda[i],")")),type="n",xaxt="n",xlab="",log="y")
axis(1,at=1:7,c(expression(g[max]),expression(g[2]),expression(g[3]),expression(g[4]),expression(g[5]),expression(g[6]),expression(g[7])),las=1)

temp_int <- HPDinterval(mcmc(sampled_variance_GA1_null[,2:8]),prob=.95)
arrows((1:7)-.15,temp_int[,1],(1:7)-.15,temp_int[,2],length=.05,code=3,angle=90,col="orange",cex=1.2)
points(c(1:7)-.15,eigen(VCV_mat_NaCl[[2]]$G1_mat/2)$values ,pch=21,bg=col_GA[1])

temp_int <- HPDinterval(mcmc(sampled_variance_GA2_null[,2:8]),prob=.95)
arrows((1:7),temp_int[,1],(1:7),temp_int[,2],length=.05,code=3,angle=90,col="orange",cex=1.2)
points(c(1:7),eigen(VCV_mat_NaCl[[3]]$G1_mat/2)$values ,pch=21,bg=col_GA[2])

temp_int <- HPDinterval(mcmc(sampled_variance_GA1_null[,2:8]),prob=.95)
arrows((1:7)+.15,temp_int[,1],(1:7)+.15,temp_int[,2],length=.05,code=3,angle=90,col="orange",cex=1.2)
points(c(1:7)+.15,eigen(VCV_mat_NaCl[[4]]$G1_mat/2)$values ,pch=21,bg=col_GA[3])

legend(3,.9,"Posterior mode",bty="n")
legend(3,.6,c("GA150","GA250","GA450"),bty="n",ncol=2,pch=21,pt.bg=col_GA,col="black")
legend(3,.3,lwd=2,"95% CI of null posterior\nmodes",bty="n",col="orange")

dev.off()


pdf(file='plots/Figure7_figure_supplement3.pdf')

plot(eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values~c(1:7),pch=16,xlim=c(1,7.5),bty="n",ylim=c(0.001,.83),ylab='Genetic variance',type="n",xaxt="n",xlab="",log="y")
axis(1,at=1:7,c(expression(g[max]),expression(g[2]),expression(g[3]),expression(g[4]),expression(g[5]),expression(g[6]),expression(g[7])),las=1)

temp_int <- HPDinterval(mcmc(sampled_variance_GA1_null_v2),prob=.95)
arrows((1:7)-.15,temp_int[,1],(1:7)-.15,temp_int[,2],length=.05,code=3,angle=90,col="orange",cex=1.2)
points(c(1:7)-.15,true_GA1_v2 ,pch=21,bg=col_GA[1])

temp_int <- HPDinterval(mcmc(sampled_variance_GA2_null_v2),prob=.95)
arrows((1:7),temp_int[,1],(1:7),temp_int[,2],length=.05,code=3,angle=90,col="orange",cex=1.2)
points(c(1:7),true_GA2_v2 ,pch=21,bg=col_GA[2])

temp_int <- HPDinterval(mcmc(sampled_variance_GA1_null_v2),prob=.95)
arrows((1:7)+.15,temp_int[,1],(1:7)+.15,temp_int[,2],length=.05,code=3,angle=90,col="orange",cex=1.2)
points(c(1:7)+.15,true_GA4_v2 ,pch=21,bg=col_GA[3])

points( (c(1:7)+0.2),eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values ,pch=21,bg="grey")

legend(3,.9,"Posterior mode",bty="n")
legend(3,.6,c("GA150","GA250","GA450","A6140"),bty="n",ncol=2,pch=21,pt.bg=c(col_GA,'gray'),col="black")
legend(3,.3,lwd=2,"95% CI of null posterior\nmodes",bty="n",col="orange")

dev.off()


