rm(list=ls());gc()
library(MCMCglmm)

# We load the data used for the G calculation
temp.env=new.env()
load('output_files/RData/Phenotypic_data.RData',envir=temp.env)
data_for_G_NaCl <- subset(temp.env$final_merged,env_label=='NaCl')
rm(temp.env);gc()

vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")
nb_trait=7

### We will run a set of MCMCglmm models with various priors - only A6140 in High Salt

i="A6140"
temp_final = subset(data_for_G_NaCl,   population == i)
### Normalize the env. factors
for(j in c('temperature',"rel_humidity","logD")){
  temp_final[,j] <- (temp_final[,j]-mean(temp_final[,j]))/sd(temp_final[,j])
}

phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))

prior.init <- list(G = list(G1 = list(V = phen.var/3, nu = 7), G2 = list(V = phen.var/3, n = 7)), 
                  R = list(V = phen.var/3, n = 7))


# flat improper prior, equivalent to REML fitting:
    prior.reml <- list(G=list(G1 = list(V = diag(7), nu = 0.00001), G2 = list(V = diag(7), nu = 0.00001)),R = list(V =diag(7), nu = 0.00001))    
# Inverse Wishart #1
    prior.iw <- list(G=list(G1 = list(V = diag(7), nu = 7), G2 = list(V = diag(7), nu = 7)),R = list(V =diag(7), nu = 7)) 
    
# Inverse Wishart  #2 - expanded
    prior.iw_exp <- list(G=list(G1 = list(V = diag(7), nu = 7 ,alpha.mu=rep(0,7),alpha.V=diag(7)*1000), G2 = list(V = diag(7), nu = 7 ,alpha.mu=rep(0,7),alpha.V=diag(7)*1000)),
                         R=list(V=diag(7),fix=1))
# Inverse Gamma
    prior.ig <- list(G=list(G1 = list(V = diag(7), nu = 6.002), G2 = list(V = diag(7), nu = 6.002)),R = list(V =diag(7), nu = 6.002))
    
    
mod.init <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32,area.F)) ~  (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                       rcov = ~us(trait):units, 
                       family = rep("gaussian", nb_trait), data = temp_final, prior = prior.init, verbose = TRUE,nitt=1500000, burnin=500000,thin=100)


mod.reml <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32,area.F)) ~  (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
rcov = ~us(trait):units, 
family = rep("gaussian", nb_trait), data = temp_final, prior = prior.reml, verbose = TRUE,nitt=1500000, burnin=500000,thin=100)

mod.iw <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32,area.F)) ~  (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                     rcov = ~us(trait):units, 
                     family = rep("gaussian", nb_trait), data = temp_final, prior = prior.iw, verbose = TRUE,nitt=1500000, burnin=500000,thin=100)

mod.iw_exp <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32,area.F)) ~  (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                     rcov = ~us(trait):units, 
                     family = rep("gaussian", nb_trait), data = temp_final, prior = prior.iw_exp, verbose = TRUE,nitt=1500000, burnin=500000,thin=100)

mod.ig <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32,area.F)) ~  (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                     rcov = ~us(trait):units, 
                     family = rep("gaussian", nb_trait), data = temp_final, prior = prior.ig, verbose = TRUE,nitt=1500000, burnin=500000,thin=100)


matrix.init = list(G1_mat = matrix(posterior.mode(mod.init$VCV)[1:nb_trait^2], nb_trait, nb_trait),VCV_Mat = mod.init$VCV)
matrix.reml = list(G1_mat = matrix(posterior.mode(mod.reml2$VCV)[1:nb_trait^2], nb_trait, nb_trait),VCV_Mat = mod.reml2$VCV)
matrix.iw = list(G1_mat = matrix(posterior.mode(mod.iw$VCV)[1:nb_trait^2], nb_trait, nb_trait),VCV_Mat = mod.iw$VCV)
matrix.iw_exp = list(G1_mat = matrix(posterior.mode(mod.iw_exp$VCV)[1:nb_trait^2], nb_trait, nb_trait),VCV_Mat = mod.iw_exp$VCV)
matrix.ig = list(G1_mat = matrix(posterior.mode(mod.ig$VCV)[1:nb_trait^2], nb_trait, nb_trait),VCV_Mat = mod.ig$VCV)

save(list=c("matrix.init","matrix.reml","matrix.iw","matrix.iw_exp","matrix.ig"),file='output_files/RData/Priors.RData')


rm(list=ls());gc()
load('output_files/RData/Priors.RData')
all_outputs = list(matrix.init,matrix.iw,matrix.reml,matrix.iw_exp,matrix.ig)
##### Plot all of them together
pdf("plots/Figure2_figure_supplement1.pdf",h=10)
par(mar=c(5,7,4,2))
vect_Var <- c(2:7,10:14,18:21,26:28,34,35,42,1,9,17,25,33,41,49)
vProb <- .95
  

vectX.init <- c(matrix.init$G1_mat/2)[vect_Var]

plot(vectX.init,c(31:11,7:1),yaxt="n",bty="n",xlim=c(-.25,.50),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)

lines(c(0,0),c(32.5,9.5))
lines(c(0,0),c(.5,6.5),col="black")
  
axis(side=1,pos=0)
  
axis(side=2,at=c(31:11,7:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF","SF*Size",
                                     "SB*FS","SB*FB","SB*BS","SB*BF","SB*Size",
                                     "FS*FB","FS*BS","FS*BF","FS*Size",
                                     "FB*BS","FB*BF","FB*Size","BS*BF","BS*Size","BF*Size",
                                     "SF","SB","FS","FB","BS","BF","Size"),las=1)
k=0
p=0
vcol=rainbow(5)
for(temp.mcmc in all_outputs){
  p=p+1
temp_95 <- HPDinterval(temp.mcmc$VCV_Mat[,1:49]/2,prob=.95)
arrows(temp_95[vect_Var,1],c(31:11,7:1)+k,temp_95[vect_Var,2],c(31:11,7:1)+k,code=3,length=.02,angle=90)

temp_80 <- HPDinterval(temp.mcmc$VCV_Mat[,1:49]/2,prob=.8)
arrows(temp_80[vect_Var,1],c(31:11,7:1)+k, temp_80[vect_Var,2],c(31:11,7:1)+k,code=3,length=0,angle=90,lwd=2,col=vcol[p])

vectX <- c(temp.mcmc$G1_mat/2)[vect_Var]
points(vectX,c(31:11,7:1)+k,pch=21,bg=vcol[p],cex=.8)

k=k-.15
}

legend(.2,25,c("Model 1","Model 2","Model 3","Model 4","Model 5"),col=vcol,lwd=2,pch=16,bty="n")
dev.off()
