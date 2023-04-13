rm(list=ls())
library(MCMCglmm)

### Here we will produce the figure 8AB

load("output_files/RData/G_matrices_high_salt.RData")
vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")

load('output_files/RData/SSCP_eigenvectors.RData')

# Load the randomized G matrices
null_matrix=new.env()
load("output_files_elife/RData_Random/Random_G_Analysis_Cemee_Pop_WI_A6140_NaCl_LIGHT.RData",envir=null_matrix)
A6140_NaCl_NULL = null_matrix$df_G1
rm(null_matrix);gc()

angle_theta <- function(x, y) {
  dot.prod <- x %*% y
  norm.x <- norm(x, type = "2")
  norm.y <- norm(y, type = "2")
  theta <- 180/pi * as.numeric(acos(dot.prod/(norm.x * norm.y)))
  as.numeric(theta)
}

gmax_NaCl <- eigen(VCV_mat_NaCl[[1]]$G1_mat)$vectors[,1]
gmax_proj_NaCl<- (t(EV_div[,1])%*%(VCV_mat_NaCl[[1]]$G1_mat/2)%*%EV_div[,1])/sum(EV_div[,1]^2)

vect_rand_piE_NaCl <- NULL
vect_rand_piE_NaCl_NULL <- NULL

vect_rand_thetaE_NaCl <- NULL
vect_rand_thetaE_NaCl_2 <- NULL
vect_rand_thetaE_NaCl_3 <- NULL

n_iter=1000
spld_idx <- sample(1:nrow(VCV_mat_NaCl[[1]]$VCV_Mat),n_iter)
k=0
for(i in spld_idx){
  k=k+1
  
  temp_NaCl <- matrix(VCV_mat_NaCl[[1]]$VCV_Mat[i,1:49],7,7)
  temp_NaCl_NULL <- matrix(A6140_NaCl_NULL[k,],7,7)
  
  temp_gmax_proj <- (t(EV_div[,1])%*%(temp_NaCl/2)%*%EV_div[,1])/sum(EV_div[,1]^2)
  temp_gmax_proj_NULL <- (t(EV_div[,1])%*%(temp_NaCl_NULL/2)%*%EV_div[,1])/sum(EV_div[,1]^2)
  vect_rand_piE_NaCl <- c(vect_rand_piE_NaCl , temp_gmax_proj/eigen(temp_NaCl/2)$values[1])	
  vect_rand_piE_NaCl_NULL <- c(vect_rand_piE_NaCl_NULL , temp_gmax_proj_NULL/eigen(temp_NaCl_NULL/2)$values[1])	

  vect_rand_thetaE_NaCl <- c(vect_rand_thetaE_NaCl,angle_theta(EV_div[,1],eigen(temp_NaCl/2)$vector[,1]))
  vect_rand_thetaE_NaCl_2 <- c(vect_rand_thetaE_NaCl_2,angle_theta(EV_div[,1],eigen(temp_NaCl/2)$vector[,2]))
  vect_rand_thetaE_NaCl_3 <- c(vect_rand_thetaE_NaCl_3,angle_theta(EV_div[,1],eigen(temp_NaCl/2)$vector[,3]))
  
}

vect_rand_thetaE_NaCl[vect_rand_thetaE_NaCl>90] = 180 - vect_rand_thetaE_NaCl[vect_rand_thetaE_NaCl>90]
vect_rand_thetaE_NaCl_2[vect_rand_thetaE_NaCl_2>90] = 180 - vect_rand_thetaE_NaCl_2[vect_rand_thetaE_NaCl_2>90]
vect_rand_thetaE_NaCl_3[vect_rand_thetaE_NaCl_3>90] = 180 - vect_rand_thetaE_NaCl_3[vect_rand_thetaE_NaCl_3>90]


Angle_rand=NULL
for(i in 1:1000){
  Angle_rand=c(Angle_rand,angle_theta(runif(7,min=(-1),max=1), runif(7,min=(-1),max=1)))
}
Angle_rand[Angle_rand>90]= 180-Angle_rand[Angle_rand>90]


pdf(file='plots_elife/Figure8AB.pdf',width=6, height=4)

layout(mat=matrix(c(1,2),1,2),w=c(1,1.3))
par(mar=c(5,4,4,2))

plot(1,1,type="n",ylim=c(0,1),xlim=c(2.2,2.38),bty="n",las=1,yaxt="n",xlab="",ylab=expression(paste("Projection along ",d[max]," (",Pi,")")),xaxt="n")
axis(2,at=c(0,0.5,1))
mtext(side=3,"A",at=2.1,cex=2)
axis(side=1,at=c(2.25,2.3),labels=c("Null","Observed"),las=2)
temp_95 <- HPDinterval(as.mcmc(vect_rand_piE_NaCl))
temp_80 <- HPDinterval(as.mcmc(vect_rand_piE_NaCl),prob=0.83)

arrows(2.3,temp_95[1,1],2.3,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(2.3,temp_80[1,1],2.3,temp_80[1,2],code=3,angle=90,length=0,col="gray",lwd=2)
points(2.3,mean(vect_rand_piE_NaCl),pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_piE_NaCl_NULL))
arrows(2.25,temp_95[1,1],2.25,temp_95[1,2],code=3,angle=90,length=0.05,col="orange",lwd=2)

par(mar=c(5,4,4,2))
plot(1,1,type="n",ylim=c(0,90),xlab="",xlim=c(1,2.5),bty="n",las=1,yaxt="n",ylab=expression(paste("Angle with ",d[max]," (",Theta,")")),xaxt="n")
axis(2,at=c(0,45,90))
mtext(side=3,"B",at=0,cex=2)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_NaCl))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_NaCl),prob=0.83)

arrows(1.5,temp_95[1,1],1.5,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(1.5,temp_80[1,1],1.5,temp_80[1,2],code=3,angle=90,length=0,col="gray",lwd=2)
points(1.5,mean(vect_rand_thetaE_NaCl),pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_NaCl_2))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_NaCl_2),prob=0.83)

arrows(1.75,temp_95[1,1],1.75,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(1.75,temp_80[1,1],1.75,temp_80[1,2],code=3,angle=90,length=0,col="gray",lwd=2)
points(1.75,mean(vect_rand_thetaE_NaCl_2),pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_NaCl_3))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_NaCl_3),prob=0.83)

arrows(2,temp_95[1,1],2,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(2,temp_80[1,1],2,temp_80[1,2],code=3,angle=90,length=0,col="gray",lwd=2)
points(2,mean(vect_rand_thetaE_NaCl_3),pch=16)

temp_95 <- HPDinterval(as.mcmc(Angle_rand))
temp_80 <- HPDinterval(as.mcmc(Angle_rand),prob=0.83)

arrows(1.1,temp_95[1,1],1.1,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(1.1,temp_80[1,1],1.1,temp_80[1,2],code=3,angle=90,length=0,col="gray",lwd=2)

axis(side=1,at=c(2,1.75,1.5,1.1),c(expression(g[3]),expression(g[2]),expression(g[max]),"Null"),las=2)

dev.off()


