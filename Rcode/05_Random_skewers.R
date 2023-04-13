rm(list = ls());gc()
library(MCMCglmm)
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
#library(parallel)
#library(RColorBrewer)


### Code from Auguirre et al. supplementary file - function to run the random skewers

#START
R.proj <- function(Gs,p,vec){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]]
  rand.vec <-matrix(,vec,n)
  for (i in 1:vec){
    b <- runif(n,-1,1)
    rand.vec[i,] <- b/(sqrt(sum(b^2)))
  }
  #generate unit length random vectors  
  proj<- function(G,b) t(b) %*% G %*% (b)
  #internal function to do projection
  G.proj <- array(,c(MCMCsamp, m, vec))
  colnames(G.proj) <- dimnames(Gs)[[3]]
  for (i in 1:vec){
    G.proj[,,i]<- t(apply(Gs, 3:4, proj, b = rand.vec[i,]))
  }
  #project each random vector through each MCMC sample of each G
  prs <- cbind(rep(1:m, each = m), 1:m) 
  prs.comp <- prs[prs[,1] < prs[,2], , drop = FALSE] 
  #setting up an index for HPD comparisons
  proj.score <-matrix(,vec,((m^2 - m)/2))
  for (k in 1:vec){
    HPD.int <- HPDinterval(as.mcmc(G.proj[,,k]), prob = p)
    proj.score[k,] <- ifelse(HPD.int[prs.comp[,1],1] > HPD.int[prs.comp[,2],2] | HPD.int[prs.comp[,2],1] > HPD.int[prs.comp[,1],2],1,0) 
  }
  #for a given random vector, examine if the HPD intervals of any pair of G matrices overlap
  vec.score <-cbind(rand.vec, proj.score)
  colnames(vec.score) <- c(1:n, paste(dimnames(Gs)[[3]][prs.comp[, 1]], ".vs.", dimnames(Gs)[[3]][prs.comp[, 2]], sep = ""))
  #collate the random vectors and the outcome of their projection on the G matrices
  sig.vec <- subset(vec.score, rowSums(vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0) 
  #collate just the random vectors that resulted in significant differences in variance
  if(dim(sig.vec)[1] <= 1) {warning("There were <= 1 significant vectors, try a larger vec or lower p"); eig.R <- "Na"}
  else{
    eig.R <- eigen(cov(sig.vec[,1:n]))
    rownames(eig.R$vectors) <- dimnames(Gs)[[1]]
    colnames(eig.R$vectors) <- c(paste("e", 1:n, sep = ""))
  }  
  #eigen analysis of the R matrix
  list(G.proj = G.proj, vec.score = vec.score, eig.R = eig.R)
}
#END


# We can load the data

load('output_files/RData/G_matrices_high_salt.RData')
vect_Pops=paste0(c("A6140","GA150","GA250","GA450"),"_NaCl")
vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")

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
#Parray <- array(, c(n, n, m, MCMCsamp))
#dimnames(Parray) <- list(traitnames, traitnames, Gnames)

#Earray1 <- array(, c(n, n, m, MCMCsamp))
#dimnames(Earray1) <- list(traitnames, traitnames, Gnames)
#Earray2 <- array(, c(n, n, m, MCMCsamp))
#dimnames(Earray2) <- list(traitnames, traitnames, Gnames)

for (i in 1:m) {
  for (j in 1:MCMCsamp) {
    Garray[, , i, j] <- matrix(MCMCarray[j, 1:(n^2), i], ncol = n)
  #  R1 <- matrix(MCMCarray[j, ((n^2) + 1):((n^2) * 2), i], ncol = n)
  #  R2 <- matrix(MCMCarray[j, (((n^2) * 2) + 1):((n^2) * 3), i], ncol = n)
  #  Earray1[, , i, j] <- R1
  #  Earray2[, , i, j] <- R2	
  #  Parray[, , i, j] <- G + R1 + R2
  }
}


# Run the random skewers analysis
MCMC.R.proj <- R.proj(Garray, p = 0.95, vec = 1000)
# Number of vectors with differences between G-matrices
table(rowSums(MCMC.R.proj$vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0 )
colSums((MCMC.R.proj$vec.score[,(n+1):(n+((m^2 - m)/2))]) ) 

# Extract the R matrix
R.structure <- lapply(MCMC.R.proj$eig.R, round, digits = 3)

#Function to do projection
proj<- function(G, b) t(b) %*% G %*% (b)


R.vec.proj <- array(, c(MCMCsamp, m, n))
for (i in 1:n){
  R.vec.proj[,,i] <- t(apply(Garray, 3:4, proj, b = MCMC.R.proj$eig.R$vectors[,i]))
}

#Genetic variance in each population in the direction of the eigenvectors of R

HPD.R.vec.proj <- array(, c(m, 2, n))
for (i in 1:n){
  HPD.R.vec.proj[,,i] <- HPDinterval(as.mcmc(R.vec.proj[,,i]), prob = 0.95)    
}
#HPD intervals for the genetic variance in each population in the direction of the eigenvectors of R
HPD.R.vec.proj.83 <- array(, c(m, 2, n))
for (i in 1:n){
  HPD.R.vec.proj.83[,,i] <- HPDinterval(as.mcmc(R.vec.proj[,,i]), prob = 0.83)    
}

# Posterior distribution for G matrices in the eigendimensions of R
R.postmode=NULL
for(i in 1:7){
  R.postmode=rbind(R.postmode,apply(R.vec.proj[,,i],2,function(x){as.numeric(posterior.mode(as.mcmc(x)))}))
}

col_GA <- c("cadetblue1", "cornflowerblue", "slateblue2")

pdf("plots/Random_Skewers.pdf",w=9)
par(mfrow=c(2,4))
for(i in 1:7){
plot(c(1:4),c(0,0,0,2),type="n",las=1,bty="n",main=paste0("e_",i),ylab="Genetic variance",xaxt="n",xlab="")
axis(side=1,at=1:4,c("A6140","GA150","GA250","GA450"))
arrows(1:4,HPD.R.vec.proj[,1,i],1:4,HPD.R.vec.proj[,2,i],code=3,angle=90,length=0.05)
arrows(1:4,HPD.R.vec.proj.83[,1,i],1:4,HPD.R.vec.proj.83[,2,i],length=0,col=c("grey",col_GA),lwd=2)
points(1:4,R.postmode[i,],pch=16)
}
dev.off()
save(list=c("R.structure","HPD.R.vec.proj","HPD.R.vec.proj.83","R.postmode"),file='output_files/RData/Random_skewers.RData')
