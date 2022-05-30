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

load('output_files/RData/Analysis_Cemee_Pop_WI_NaCl.RData')

true_G_NGM = VCV_mat_NGM[[1]]/2
true_G_NaCl = VCV_mat_NaCl[[1]]/2


plot_G_coefs <- function(true_G_NGM, true_G_NaCl,pop_name){
  
  par(mar=c(5,7,4,2))
  vect_Var <- c(2:6,9:12,16:18,23,24,30,1,8,15,22,29,36)
  vProb <- .95
  
  vectX_NGM <- c(true_G_NGM$G1_mat/2)[vect_Var]
  vectX_NaCl <- c(true_G_NaCl$G1_mat/2)[vect_Var]
  
  plot(vectX_NGM,c(24:10,6:1),yaxt="n",bty="n",xlim=c(-.25,.50),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)
  #mtext(side=2,"Phenotypic traits    \n Diagonal                                         Off-diagonal                     ",padj=-2,cex=1.2)
  lines(c(0,0),c(24.5,8.5))
  lines(c(0,0),c(.5,5.5),col="red")
  
  axis(side=1,pos=0)
  
  axis(side=2,at=c(24:10,6:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                                       "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                                       "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                                       "SF","SB","FS","FB","BS","BF"),las=1)
  
  temp_95 <- HPDinterval(true_G_NGM$VCV_Mat[,1:36]/2,prob=.95)
  arrows(temp_95[vect_Var,1],c(24:10,6:1),temp_95[vect_Var,2],c(24:10,6:1),code=3,length=.02,angle=90)
  
  temp_95 <- HPDinterval(true_G_NaCl$VCV_Mat[,1:36]/2,prob=.95)
  arrows(temp_95[vect_Var,1],c(24:10,6:1)+.3,temp_95[vect_Var,2],c(24:10,6:1)+.3,code=3,length=.02,angle=90)
  
  temp_80 <- HPDinterval(true_G_NGM$VCV_Mat[,1:36]/2,prob=.8)
  arrows(temp_80[vect_Var,1],c(24:10,6:1), temp_80[vect_Var,2],c(24:10,6:1),code=3,length=0,angle=90,lwd=2,col="cornflowerblue")
  
  temp_80 <- HPDinterval(true_G_NaCl$VCV_Mat[,1:36]/2,prob=.8)
  arrows(temp_80[vect_Var,1],c(24:10,6:1)+.3, temp_80[vect_Var,2],c(24:10,6:1)+.3,code=3,length=0,angle=90,lwd=2,col="firebrick3")
  
  points(vectX_NGM,c(24:10,6:1),pch=21,bg="black",cex=.6)
  points(vectX_NaCl,c(24:10,6:1)+.3,pch=21,bg="black",cex=.6)
  
  
  legend(.1,23,c("Low Salt","High Salt"),lwd=2,bty="n",col=c("cornflowerblue","firebrick3"))
  
}

pdf(file='plots/G_mat_NGM_NaCl_A6140.pdf',h=8,w=5.5)
plot_G_coefs(VCV_mat_NGM[[1]],VCV_mat_NaCl[[1]],"A6140")
dev.off()

######################################################################
##### Angle between the G-max of the NGM/NaCl matrices     ####
######################################################################

angle_theta <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- 180/pi * as.numeric(acos(dot.prod/(norm.x * norm.y)))
	if(theta>90) theta <- (180-theta)
	as.numeric(theta)
}

vect_NGM_sampled <- sample(1:nrow(VCV_mat_NGM[[1]]$VCV_Mat/2),1000)
vect_NaCl_sampled <- sample(1:nrow(VCV_mat_NaCl[[1]]$VCV_Mat/2),1000)
all_angles <- NULL
for(i in 1:1000){
	temp_NGM_mat <- matrix(VCV_mat_NGM[[1]]$VCV_Mat[vect_NGM_sampled[i],1:36],6,6)/2
	temp_NaCl_mat <- matrix(VCV_mat_NaCl[[1]]$VCV_Mat[vect_NaCl_sampled[i],1:36],6,6)/2	
	all_angles <- c(all_angles ,angle_theta(eigen(temp_NGM_mat)$vectors[,1],eigen(temp_NaCl_mat)$vectors[,1]))
	
}

angle_theta(eigen(VCV_mat_NGM[[1]]$G1_mat)$vectors[,1],eigen(VCV_mat_NaCl[[1]]$G1_mat)$vectors[,1])
HPDinterval(as.mcmc(all_angles)) 

######################################################################
##### Eigen decomposition of the NGM/NaCl matrices     ####
######################################################################

sampled_variance_NGM=NULL
sampled_variance_NaCl=NULL

sampled_variance_GA1=NULL
sampled_variance_GA2=NULL
sampled_variance_GA4=NULL

sampled_variance_GA1_NGM=NULL
sampled_variance_GA2_NGM=NULL
sampled_variance_GA4_NGM=NULL

for(i in 1:nrow(VCV_mat_NGM[[1]]$VCV_Mat)){
  sampled_variance_NGM =rbind(sampled_variance_NGM,c(sum(eigen(matrix(VCV_mat_NGM[[1]]$VCV_Mat[i,1:36],6,6)/2)$values),eigen(matrix(VCV_mat_NGM[[1]]$VCV_Mat[i,1:36],6,6)/2)$values))
  sampled_variance_NaCl=rbind(sampled_variance_NaCl,c(sum(eigen(matrix(VCV_mat_NaCl[[1]]$VCV_Mat[i,1:36],6,6)/2)$values),eigen(matrix(VCV_mat_NaCl[[1]]$VCV_Mat[i,1:36],6,6)/2)$values))
  
  sampled_variance_GA1 =rbind(sampled_variance_GA1,c(sum(eigen(matrix(VCV_mat_NaCl[[2]]$VCV_Mat[i,1:36],6,6)/2)$values),eigen(matrix(VCV_mat_NaCl[[2]]$VCV_Mat[i,1:36],6,6)/2)$values))
  sampled_variance_GA2 =rbind(sampled_variance_GA2,c(sum(eigen(matrix(VCV_mat_NaCl[[3]]$VCV_Mat[i,1:36],6,6)/2)$values),eigen(matrix(VCV_mat_NaCl[[3]]$VCV_Mat[i,1:36],6,6)/2)$values))
  sampled_variance_GA4 =rbind(sampled_variance_GA4,c(sum(eigen(matrix(VCV_mat_NaCl[[4]]$VCV_Mat[i,1:36],6,6)/2)$values),eigen(matrix(VCV_mat_NaCl[[4]]$VCV_Mat[i,1:36],6,6)/2)$values))
  
  sampled_variance_GA1_NGM =rbind(sampled_variance_GA1_NGM,c(sum(eigen(matrix(VCV_mat_NGM[[2]]$VCV_Mat[i,1:36],6,6)/2)$values),eigen(matrix(VCV_mat_NGM[[2]]$VCV_Mat[i,1:36],6,6)/2)$values))
  sampled_variance_GA2_NGM =rbind(sampled_variance_GA2_NGM,c(sum(eigen(matrix(VCV_mat_NGM[[3]]$VCV_Mat[i,1:36],6,6)/2)$values),eigen(matrix(VCV_mat_NGM[[3]]$VCV_Mat[i,1:36],6,6)/2)$values))
  sampled_variance_GA4_NGM =rbind(sampled_variance_GA4_NGM,c(sum(eigen(matrix(VCV_mat_NGM[[4]]$VCV_Mat[i,1:36],6,6)/2)$values),eigen(matrix(VCV_mat_NGM[[4]]$VCV_Mat[i,1:36],6,6)/2)$values))
  
}



pdf(file='plots/G_A6140_Variance_decomposition.pdf',h=9,w=4)

layout(mat=matrix(c(1,2),nrow=2), widths = c(1,1),heights=c(5,6))


plot(c(0,0.2),c(
  sum(eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values),
  sum(eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values)),bg=c("black","black"),pch=21,ylab='Total genetic variance',xlim=c(-.2,.4),ylim=c(0,.8),type="n",bty="n",xaxt="n",xlab="")

axis(1,at=c(0,.2),c("Low salt","High salt"),las=2)
temp_int <- HPDinterval(mcmc(sampled_variance_NGM[,1]),prob=.95)
arrows(0,temp_int[1],0,temp_int[2],length=.05,code=3,angle=90)
temp_int <- HPDinterval(mcmc(sampled_variance_NaCl[,1]),prob=.95)
arrows(0.2,temp_int[1],0.2,temp_int[2],length=.05,code=3,angle=90)

temp_int <- HPDinterval(mcmc(sampled_variance_NGM[,1]),prob=.8)
arrows(0,temp_int[1],0,temp_int[2],length=0,code=3,angle=90,lwd=2,col="cornflowerblue")
temp_int <- HPDinterval(mcmc(sampled_variance_NaCl[,1]),prob=.8)
arrows(0.2,temp_int[1],0.2,temp_int[2],length=0,code=3,angle=90,lwd=2,col="firebrick3")

points(c(0,0.2),c(
  sum(eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values),
  sum(eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values)),bg=c("black","black"),pch=21)

###########

plot(eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values~c(1:6),pch=16,xlim=c(1,6.5),bty="n",ylim=c(0,.8),ylab='Genetic variance',type="n",xaxt="n",xlab="")
axis(1,at=1:6,c(expression(lambda[max]),expression(lambda[2]),expression(lambda[3]),expression(lambda[4]),expression(lambda[5]),expression(lambda[6])),las=1)
for(i in 1:6){
  temp_int <- HPDinterval(mcmc(sampled_variance_NGM[,(i+1)]),prob=.95)
  arrows(i,temp_int[1],i,temp_int[2],length=.05,code=3,angle=90)
  temp_int <- HPDinterval(mcmc(sampled_variance_NaCl[,(i+1)]),prob=.95)
  arrows(i+0.2,temp_int[1],i+0.2,temp_int[2],length=.05,code=3,angle=90)
  
  temp_int <- HPDinterval(mcmc(sampled_variance_NGM[,(i+1)]),prob=.8)
  arrows(i,temp_int[1],i,temp_int[2],length=0,code=3,angle=90,lwd=2,col="cornflowerblue")
  temp_int <- HPDinterval(mcmc(sampled_variance_NaCl[,(i+1)]),prob=.8)
  arrows(i+0.2,temp_int[1],i+0.2,temp_int[2],length=0,code=3,angle=90,lwd=2,col="firebrick3")
  
}
points( (c(1:6)),eigen(VCV_mat_NGM[[1]]$G1_mat/2)$values ,pch=16,col="black")
points( (c(1:6)+0.2),eigen(VCV_mat_NaCl[[1]]$G1_mat/2)$values ,pch=16,col="black")

dev.off()


####
all_angles_GA1 <- NULL
all_angles_GA2 <- NULL
all_angles_GA4 <- NULL

vect_GA1_sampled <- sample(1:nrow(VCV_mat_NaCl[[2]]$VCV_Mat/2),1000)
vect_GA2_sampled <- sample(1:nrow(VCV_mat_NaCl[[3]]$VCV_Mat/2),1000)
vect_GA4_sampled <- sample(1:nrow(VCV_mat_NaCl[[4]]$VCV_Mat/2),1000)

vect_GA1_NGM_sampled <- sample(1:nrow(VCV_mat_NaCl[[2]]$VCV_Mat/2),1000)
vect_GA2_NGM_sampled <- sample(1:nrow(VCV_mat_NaCl[[3]]$VCV_Mat/2),1000)
vect_GA4_NGM_sampled <- sample(1:nrow(VCV_mat_NaCl[[4]]$VCV_Mat/2),1000)

all_angles_GA1_NGM <- NULL
all_angles_GA2_NGM <- NULL
all_angles_GA4_NGM <- NULL

for(i in 1:1000){

	temp_NaCl_mat <- matrix(VCV_mat_NaCl[[1]]$VCV_Mat[vect_NaCl_sampled[i],1:36],6,6)/2	
	temp_NGM_mat <-  matrix(VCV_mat_NGM[[1]]$VCV_Mat[vect_NGM_sampled[i],1:36],6,6)/2

	temp_GA1_mat  <- matrix(VCV_mat_NaCl[[2]]$VCV_Mat[vect_GA1_sampled[i],1:36],6,6)/2
	temp_GA2_mat  <- matrix(VCV_mat_NaCl[[3]]$VCV_Mat[vect_GA2_sampled[i],1:36],6,6)/2
	temp_GA4_mat  <- matrix(VCV_mat_NaCl[[4]]$VCV_Mat[vect_GA4_sampled[i],1:36],6,6)/2
	
	all_angles_GA1 <- c(all_angles_GA1 ,angle_theta(eigen(temp_GA1_mat)$vectors[,1],eigen(temp_NaCl_mat)$vectors[,1]))
	all_angles_GA2 <- c(all_angles_GA2 ,angle_theta(eigen(temp_GA2_mat)$vectors[,1],eigen(temp_NaCl_mat)$vectors[,1]))	
	all_angles_GA4 <- c(all_angles_GA4 ,angle_theta(eigen(temp_GA4_mat)$vectors[,1],eigen(temp_NaCl_mat)$vectors[,1]))

	temp_GA1_NGM_mat  <- matrix(VCV_mat_NGM[[2]]$VCV_Mat[vect_GA1_NGM_sampled[i],1:36],6,6)/2
	temp_GA2_NGM_mat  <- matrix(VCV_mat_NGM[[3]]$VCV_Mat[vect_GA2_NGM_sampled[i],1:36],6,6)	/2
	temp_GA4_NGM_mat  <- matrix(VCV_mat_NGM[[4]]$VCV_Mat[vect_GA4_NGM_sampled[i],1:36],6,6)/2
	
	all_angles_GA1_NGM <- c(all_angles_GA1_NGM ,angle_theta(eigen(temp_GA1_NGM_mat)$vectors[,1],eigen(temp_NGM_mat)$vectors[,1]))
	all_angles_GA2_NGM <- c(all_angles_GA2_NGM ,angle_theta(eigen(temp_GA2_NGM_mat)$vectors[,1],eigen(temp_NGM_mat)$vectors[,1]))	
	all_angles_GA4_NGM <- c(all_angles_GA4_NGM ,angle_theta(eigen(temp_GA4_NGM_mat)$vectors[,1],eigen(temp_NGM_mat)$vectors[,1]))
	
	}

save(list=ls(),file='output_files/RData/Plot_all_G_matrices_no_PLAST.RData')
rm(list=ls())
######################################################################
######################################################################


v_col = c("chartreuse","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")

pdf(file='plots/G_mat_NaCl_Evolution.pdf',h=8,w=5.5)


par(mar=c(5,7,4,2))
vect_Var <- c(2:6,9:12,16:18,23,24,30,1,8,15,22,29,36)
vProb <- .95

vectX_NGM <- c(true_G_NGM$G1_mat)[vect_Var]
vectX_NaCl <- c(true_G_NaCl$G1_mat)[vect_Var]

plot(c(VCV_mat_NaCl[[1]]$G1_mat/2)[vect_Var],c(24:10,6:1),yaxt="n",bty="n",xlim=c(-.25,.50),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)

lines(c(0,0),c(24.5,8.5))
lines(c(0,0),c(.5,5.5),col="red")

axis(side=1,pos=0)

axis(side=2,at=c(24:10,6:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
"SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
"FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
"SF","SB","FS","FB","BS","BF"),las=1)

for(i in 1:4){

temp_95 <- HPDinterval(VCV_mat_NaCl[[i]]$VCV_Mat[,1:36]/2,prob=.95)
arrows(temp_95[vect_Var,1],c(24:10,6:1)+(.15*(i-1)),temp_95[vect_Var,2],c(24:10,6:1)+(.15*(i-1)),code=3,length=.02,angle=90)

temp_80 <- HPDinterval(VCV_mat_NaCl[[i]]$VCV_Mat[,1:36]/2,prob=.8)
arrows(temp_80[vect_Var,1],c(24:10,6:1)+(.15*(i-1)),
temp_80[vect_Var,2],c(24:10,6:1)+(.15*(i-1)),code=3,length=0,angle=90,lwd=2,col=c("orange",v_col[2:4])[i])

points(c(VCV_mat_NaCl[[i]]$G1_mat/2)[vect_Var],c(24:10,6:1)+(.15*(i-1)),pch=21,bg="black",cex=.6)
}

legend(.15,23,c("A6140","GA150","GA250","GA450"),lwd=2,bty="s",col=c("orange",v_col[2:4]))

dev.off()


save(list=ls(),file='output_files/RData/Plot_all_G_matrices_no_PLAST.RData')






