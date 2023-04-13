rm(list = ls())
gc()
library(MCMCglmm)
library(parallel)

#library(psych)
#library(ggplot2)
#library(dplyr)
#library(gplots)
#library(data.table)
#library(matrixStats)
#library(boot)
#library(Rmisc)
#library(dae)
#library(nlme)

#library(RColorBrewer)

angle_eigenV <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- acos(dot.prod/(norm.x * norm.y))
	as.numeric(theta)
}


load('output_files/RData/G_matrices_high_salt.RData')
vect_Pops=paste0(c("A6140","GA150","GA250","GA450"),"_NaCl")
vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")

####
MCMCtot <- nrow(VCV_mat_NaCl[[1]]$VCV_Mat)
MCMCsamp <- 1000 
n <- 7 #number of traits
m <- 4 #number of matrices to compare
r <- 3 #number of random effects specified in the model.
traitnames <- vect_P_traits #trait names
Gnames <- vect_Pops

MCMCarray <- array(, c(MCMCsamp, (n^2) * r, m)) #empty array
MCMCarray[, , 1] <- as.matrix(VCV_mat_NaCl[[1]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 2] <- as.matrix(VCV_mat_NaCl[[2]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 3] <- as.matrix(VCV_mat_NaCl[[3]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 4] <- as.matrix(VCV_mat_NaCl[[4]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])

Garray <- array(, c(n, n, m, MCMCsamp))
dimnames(Garray) <- list(traitnames, traitnames, Gnames)
Parray <- array(, c(n, n, m, MCMCsamp))
dimnames(Parray) <- list(traitnames, traitnames, Gnames)

Earray1 <- array(, c(n, n, m, MCMCsamp))
dimnames(Earray1) <- list(traitnames, traitnames, Gnames)
Earray2 <- array(, c(n, n, m, MCMCsamp))
dimnames(Earray2) <- list(traitnames, traitnames, Gnames)

for (i in 1:m) {
	for (j in 1:MCMCsamp) {
		G <- matrix(MCMCarray[j, 1:(n^2), i], ncol = n)
		R1 <- matrix(MCMCarray[j, ((n^2) + 1):((n^2) * 2), i], ncol = n)
		R2 <- matrix(MCMCarray[j, (((n^2) * 2) + 1):((n^2) * 3), i], ncol = n)
		Garray[, , i, j] <- G
		Earray1[, , i, j] <- R1
		Earray2[, , i, j] <- R2	
		Parray[, , i, j] <- G + R1 + R2
	}
}

source('Rcode/Appendix_Tensor/02_functions_tensor.R', chdir = TRUE)

HHGarray <- array(, c(n, n, m, MCMCsamp))
for (k in 1:MCMCsamp) {
	for (j in 1:m) {
		P <- inv.rootP(Parray[, , j, k])
		HHGarray[, , j, k] <- P %*% Garray[, , j, k] %*% P
	}
}

load("output_files/RData/Phenotypic_data.RData")
data_for_G_NaCl2= subset(final_merged,env_label=="NaCl" & population %in%c("A6140",paste0("GA",c(1,2,4),"50")))
data_for_G_NaCl2$population <- paste0(data_for_G_NaCl2$population,"_NaCl")
data_for_G_NaCl2$pop_label <- paste0(data_for_G_NaCl2$pop_label,"_NaCl")

## Create a pedigree
df_for_tensor= data_for_G_NaCl2
length(unique(df_for_tensor$pop_label)) # 350

ped_all = rbind(
#first all the lines
data.frame(id = as.character(unique(df_for_tensor$pop_label)), dam = NA, sire = NA,stringsAsFactors=FALSE),
#then all the phenotyped lines
data.frame(id = 1:nrow(df_for_tensor), dam = as.character(df_for_tensor$pop_label), sire = as.character(df_for_tensor$pop_label),stringsAsFactors=FALSE)
)
for(i in 1:3) ped_all[,i]=as.factor(ped_all[,i])

population_for_ped <- data.frame(pop_label= c(as.character(unique(df_for_tensor$pop_label)),as.character(df_for_tensor$pop_label)),stringsAsFactors=FALSE)

population_for_ped$population=NA
for(i in 1:nrow(population_for_ped)) population_for_ped$population[i] = as.character(subset(unique(df_for_tensor[,c("population","pop_label")]),pop_label==population_for_ped$pop_label[i])$population)

rand.Garray <- array(, c(n, n, m, MCMCsamp))
rand.Garray_corrected <- array(, c(n, n, m, MCMCsamp))
dimnames(rand.Garray) <- list(traitnames, traitnames, Gnames)

df_for_tensor$population=as.factor(as.character(df_for_tensor$population))
rm(i)

# Here we save a file that could be used on a server to compute the randomized
# eigentensors that are computationaly demanding.
save(list=ls(),file='output_files/RData/File_for_parallel_processing_Tensor.RData')

run_parallel_MCMC <- function(i){

	library(MCMCglmm)
	library(dae)
  library(data.table)
	load('output_files/RData/File_for_parallel_processing_Tensor.RData')
	rand.Garray_corrected_parallel <- array(, c(n, n, m))
	nb_trait=7

##### NaCl - Generate phenotypic values from the pedigree
	A6140_NaCl.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[1]), Garray[, , 1, i]/2)
	GA150_NaCl.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[2]), Garray[, , 2, i]/2)
	GA250_NaCl.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[3]), Garray[, , 3, i]/2)
	GA450_NaCl.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[4]), Garray[, , 4, i]/2)

	
	a.pop <- cumsum(as.numeric(tapply(ped_all$id[!is.na(ped_all$dam)], population_for_ped$population[!is.na(ped_all$dam)], length)))
	pop.bv <- rbind(A6140_NaCl.bv, GA150_NaCl.bv, GA250_NaCl.bv, GA450_NaCl.bv)
  
	#
	rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1], replace = F), ]
	
	## Here we have to compute the random Garray using the morissey technique and save them in a list
	
	ped_lines_all <- subset(ped_all,!is.na(dam))
	rand.model_MCMC=list()
	for(k in 1:4){

	k_pop=Gnames[k]
	ped_lines_current <- subset(ped_lines_all, df_for_tensor$population==k_pop)
	
	vect_sire_labels <- paste0(tstrsplit(ped_all$id,"L")[[1]]
	,"_",tstrsplit(ped_all$id,"_")[[2]])
	sire.bvs <- rand.pop.bv[vect_sire_labels ==k_pop & is.na(ped_all$dam),]
	sire.bvs <- cbind(data.frame(sire=ped_all$id[vect_sire_labels==k_pop & is.na(ped_all$dam)]),sire.bvs)
	
	ped_lines_current <- merge(ped_lines_current,sire.bvs)[,c(1,4:10)]
	#Vectors of E variance	
	z <- t(apply(ped_lines_current[,2:8],1,function(x){rmvnorm(x,Earray2[,,k,i])}))
	ped_lines_current[,2:8] <- z
	names(ped_lines_current) <- c("pop_label",traitnames)
	phen.var = diag(nb_trait) * diag(var(z))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait)), R = list(V = phen.var/3, n = nb_trait))


	rand.model_MCMC.temp <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32,area.F)) ~  trait - 1, random = ~us(trait):pop_label ,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = ped_lines_current, prior = prior_mod, verbose = FALSE,nitt=150000, burnin=50000,thin=100)
	
		rand.Garray_corrected_parallel[,,k] <- matrix(posterior.mode(rand.model_MCMC.temp$VCV)[1:49], ncol = n)
	}
	
return(rand.Garray_corrected_parallel)	
}


clust <- makeCluster(22)
param_list=list()
for(i in 1:1000) param_list[[i]] <- i

List_output <-parLapply(clust, param_list , run_parallel_MCMC)

save(list="List_output",file="output_files/RData/Tensor_processed.Rdata")

### End of the parallel thread

rm(list=ls())
gc()
#### Back on local computer
load('output_files/RData/File_for_parallel_processing_Tensor.RData')
load('output_files/RData/Tensor_processed.Rdata')
load('output_files/RData/G_matrices_high_salt.RData')

source('Rcode/Appendix_Tensor/02_functions_tensor.R', chdir = TRUE)

for(i in 1:MCMCsamp){
	for(k in 1:m){
		 rand.Garray_corrected[,,k,i] <- matrix(List_output[[i]][,,k], ncol = n)
	}
}
dimnames(rand.Garray_corrected) <- list(traitnames, traitnames, Gnames)
MCMC.covtensor <- covtensor(Garray)

nnonzero <- min(n * (n + 1)/2, m - 1)
MCMC.covtensor.rand <- covtensor(rand.Garray_corrected)

HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[, 1:nnonzero]), prob = 0.95), 
HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[, 1:nnonzero]), prob = 0.95))
round(HPD.eT.val, 3)


HPD.eT.val_80 <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[, 1:nnonzero]), prob = 0.80), 
HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[, 1:nnonzero]), prob = 0.80))
round(HPD.eT.val_80, 3)

#
# Figure A1

pdf("plots/Tensor_NaCl_A1.pdf")
par(mfrow=c(1,1))

plot((1:nnonzero)-0.2,unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),xlab="Eigentensor",ylab=expression(paste("Genetic divergence (",alpha,")")),pch=16,cex=1,xaxt="n",frame.plot=F,xlim=c(0.5,3.5),ylim=c(0,max(HPD.eT.val)),type="n")
axis(1,at=1:nnonzero,labels=c(paste("E",rep(1:nnonzero),sep="")))

arrows((1:nnonzero)-0.2, HPD.eT.val[,2],(1:nnonzero)-0.2,HPD.eT.val[,1],length=0.1,angle=90,code=3)
arrows((1:nnonzero)-0.2, HPD.eT.val_80[,2],(1:nnonzero)-0.2, HPD.eT.val_80[,1],length=0.1,angle=90,code=0,lwd=4,col="grey")

arrows((1:nnonzero)+0.2, HPD.eT.val[,4],(1:nnonzero)+0.2,HPD.eT.val[,3],length=0.1,angle=90,lty=5,code=3)

arrows((1:nnonzero)+0.2, HPD.eT.val_80[,4],(1:nnonzero)+0.2, HPD.eT.val_80[,3],length=0.1,angle=90,code=0,lwd=4,col="orange")

points((1:nnonzero)-0.2,unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),pch=21,cex=1,bg="black")

points((1:nnonzero)+0.2, unique(MCMC.covtensor.rand$tensor.summary[1:(n*nnonzero),1]),pch=8,cex=1,col="grey")

legend(1.5,.23,legend=c("observed","randomised"),lty=c(1,5),pch=c(16,8),cex=1,bty="n",col=c("grey","orange"))
dev.off()

HPD.tensor.coord <- array(,c(m,2,nnonzero))
dimnames(HPD.tensor.coord) <- list(Gnames,c("lower","upper"), paste("E",1:3,sep=" "))
for (i in 1:m){
  for (j in 1:nnonzero){
    HPD.tensor.coord[i,,j] <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[i,j,]),prob=0.95)[1:2]
  }
}


HPD.tensor.coord_80 <- array(,c(m,2,nnonzero))
dimnames(HPD.tensor.coord_80) <- list(Gnames,c("lower","upper"), paste("E",1:3,sep=" "))
for (i in 1:m){
  for (j in 1:nnonzero){
    HPD.tensor.coord_80[i,,j] <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[i,j,]),prob=0.83)[1:2]
  }
}

col_GA <- c("cadetblue1", "cornflowerblue", "slateblue2")
#Figure A2 / Only E1 here
pdf('plots/Tensor_NaCl_A2.pdf')
k=1
par(mfrow=c(1,1))
par(mar=c(5,4,4,2))

plot(1:m,MCMC.covtensor$av.G.coord[,k,],ylab="",xlab="",pch=16,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(-1.5,.1),main = "",type="n")
axis(1,at=1:m,labels=tstrsplit(Gnames,"_")[[1]],las=2)
arrows(1:m,HPD.tensor.coord[,2,k],1:m,HPD.tensor.coord[,1,k],length=0.05,angle=90,code=3)
arrows(1:m,HPD.tensor.coord_80[,2,k],1:m, HPD.tensor.coord_80[,1,k],length=0.05,angle=90,code=0,col=c("grey",col_GA),lwd=3)

points(1:m,MCMC.covtensor$av.G.coord[,k,],pch=21,bg="grey")
mtext(dimnames(MCMC.covtensor$av.G.coord)[[2]][k],side=3,at=0,font=2)

dev.off()

round(MCMC.covtensor$tensor.summary[1:(n*2),2:dim(MCMC.covtensor$tensor.summary)[2]], 3)

# How much variation is explained ?

abs(MCMC.covtensor$tensor.summary[1,2])/sum(abs(MCMC.covtensor$tensor.summary[1:n,2])) # 83 %
#e11: 90%

e11 <- c(as.numeric(MCMC.covtensor$tensor.summary[1,3:dim(MCMC.covtensor$tensor.summary)[2]]))
e11.proj <- apply(Garray, 3:4, proj, b = e11)
HPD.e11 <- HPDinterval(t(as.mcmc(e11.proj)),prob = 0.95)
HPD.e11_80 <- HPDinterval(t(as.mcmc(e11.proj)),prob = 0.80)

pdf("plots/Tensor_NaCl_A3.pdf")
par(mfrow=c(1,1))
par(mar=c(5,4,4,2))
plot(1:m,rowMeans(e11.proj),ylab=expression(paste("Genetic variance (",lambda,")")),xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e11))),las=2,type="n")
axis(1,at=1:m,labels=c("A6140","GA150","GA250","GA450"),las=2)
arrows(1:m,HPD.e11[,2],1:m,HPD.e11[,1],length=0.1,angle=90,code=3)
arrows(1:m,HPD.e11_80[,2],1:m,HPD.e11_80[,1],length=0.1,angle=90,code=0,lwd=3,col=c("grey",col_GA))
points(1:m,rowMeans(e11.proj),pch=21,bg="grey")
mtext("e11 (90%)",side=3,at=0,font=2)
dev.off()

## Angle between e11 and emax
angle_theta <- function(x, y) {
  dot.prod <- x %*% y
  norm.x <- norm(x, type = "2")
  norm.y <- norm(y, type = "2")
  theta <- 180/pi * as.numeric(acos(dot.prod/(norm.x * norm.y)))
  as.numeric(theta)
}

temp.env=new.env()
load('output_files/RData/Random_skewers.RData',envir=temp.env)
emax=temp.env$R.structure$vectors[,1]
rm(temp.env);gc()

angle_theta(emax,e11)

save(list=ls(),file='output_files/RData/Tensor_processed.Rdata')
