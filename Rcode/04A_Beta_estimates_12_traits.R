rm(list=ls())
library(MCMCglmm)
library(lme4)
library(data.table)


load('output_files/RData/Analysis_Cemee_Pop_WI_NaCl.RData')
load("data/fertility.rda")

# Load the fecundity data from a previous article
mean_relative_fit <- subset(ecoefs, pop=="A6140" & env=="NaCl")
# Compute fecundity as 'relative fitness'
mean_relative_fit$w <- exp(mean_relative_fit$fertility)/mean(exp(mean_relative_fit$fertility))

# Concatenate low and high salt measurements
data_for_G_NGM_A6140 <- subset(data_for_G_NGM,population=="A6140")
data_for_G_NaCl_A6140 <- subset(data_for_G_NaCl,population=="A6140")

names(data_for_G_NGM_A6140)[names(data_for_G_NGM_A6140)%in%vect_P_traits] <- paste0(vect_P_traits,"_NGM")
names(data_for_G_NaCl_A6140)[names(data_for_G_NaCl_A6140)%in%vect_P_traits] <- paste0(vect_P_traits,"_NaCl")

test_unique_NGM = paste(data_for_G_NGM_A6140$pop_label,data_for_G_NGM_A6140$date_str)
test_unique_NaCl = paste(data_for_G_NaCl_A6140$pop_label,data_for_G_NaCl_A6140$date_str)

data_for_G_NGM_A6140  = data_for_G_NGM_A6140[order(test_unique_NGM),]
data_for_G_NaCl_A6140 = data_for_G_NaCl_A6140[order(test_unique_NaCl),]

shared_uniq = test_unique_NGM[test_unique_NGM%in%test_unique_NaCl]

data_for_G_NGM_A6140 = subset(data_for_G_NGM_A6140, test_unique_NGM%in%test_unique_NaCl)
data_for_G_NaCl_A6140 = subset(data_for_G_NaCl_A6140, test_unique_NaCl%in%test_unique_NGM)

names(data_for_G_NGM_A6140)[names(data_for_G_NGM_A6140)%in%c("temperature","logD","rel_humidity")]=
    paste0(names(data_for_G_NGM_A6140)[names(data_for_G_NGM_A6140)%in%c("temperature","logD","rel_humidity")],"_NGM")
names(data_for_G_NaCl_A6140)[names(data_for_G_NaCl_A6140)%in%c("temperature","logD","rel_humidity")]=
  paste0(names(data_for_G_NaCl_A6140)[names(data_for_G_NaCl_A6140)%in%c("temperature","logD","rel_humidity")],"_NaCl")

env_covariates = paste0(rep(c("temperature","logD","rel_humidity"),2),rep(c("_NGM","_NaCl"),each=3))
vect_P_traits_all = paste0(rep(vect_P_traits,2),rep(c("_NGM","_NaCl"),each=6)) 

retained_names=c(vect_P_traits_all,env_covariates,'date_str',"pop_label","is_2012")

all_data <- merge(data_for_G_NGM_A6140[,names(data_for_G_NGM_A6140)%in%retained_names],
                  data_for_G_NaCl_A6140[,names(data_for_G_NaCl_A6140)%in%retained_names])

## Final data with all traits (6 TR in 2 environments)
dim(all_data) # 342

## Add the relative fitness 
mean_relative_fit <- mean_relative_fit[,c("line","w")]
names(mean_relative_fit)[1] = 'pop_label'

dim(all_data)
all_data_with_w <- (merge(all_data, mean_relative_fit))
dim(all_data_with_w) # 320

vect_P_traits_all=c(vect_P_traits_all,"w")
nb_trait = 13
k=0

	phen.var = diag(nb_trait) * diag(var(subset(all_data_with_w, select = vect_P_traits_all)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC <- MCMCglmm(cbind(c(T12_NGM, T13_NGM, T21_NGM, T23_NGM, T31_NGM, T32_NGM,T12_NaCl, T13_NaCl, T21_NaCl, T23_NaCl, T31_NaCl, T32_NaCl,w)) ~ 
	                         (temperature_NGM+rel_humidity_NGM+logD_NGM)^3 + (temperature_NaCl+rel_humidity_NaCl+logD_NaCl)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = all_data_with_w, prior = prior_mod, verbose = TRUE,nitt=500000, burnin=50000,thin=100)
	
	pdf('output_files/auto_corr_plots_MCMCglmm/GZW_all12traits.pdf',h=20,w=20)
	par(mfrow = c(7, 2), mar = c(2, 2, 1, 0))
	plot(model_MCMC$Sol, auto.layout = F,ylim=c(-1,1))
	plot.acfs(model_MCMC$Sol)
  dev.off()
  
  post_dist = posterior.mode(model_MCMC$VCV)

	VCV_with_w_NGM_12t <- list(Population = "A6_with_12t", N_measurement = nrow(temp_final), G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)

	VCV_with_w_NGM_12t$G1_mat


pdf(file='plots/GZW_all_12t.pdf',h=8,w=5.5)
	
	par(mar=c(5,7,4,2))
	vect_Var <- c(2:6,9:12,16:18,23,24,30,1,8,15,22,29,36)
	
	vect_12t <- c(2:13,16:26,30:39,44:52,58:65,72:78,
	              86:91,100:104,114:117,128:130,142:143,156,14*(1:13)-13)
	
	vect_12t_v2 <- c(c(2:6,16:19,30:32,44,45,58),c(2:6,16:19,30:32,44,45,58)+84,c(157:168),c(14*(1:13)-13))
	(15+15+12+13)
	
	vProb <- .95
	
	vectX_12t <- c(VCV_with_w_NGM_12t$G1_mat)[vect_12t]
	vectX_12t_v2 <- c(VCV_with_w_NGM_12t$G1_mat)[vect_12t_v2]
	length(vect_12t)
	
	length(c(29:34))
	vect_y_v2 =c(c(11:25),c(11:25)+.2,c(1:6),c(1:6)+.2,c(29:34),c(29:34)+.2,7)
	vect_col_v2=(c(rep("lightblue",15),rep("firebrick3",15),rep("darkgreen",12),rep("lightblue",6),rep("firebrick3",6),"darkgreen"))
	
	plot(vectX_12t_v2,vect_y_v2,ylim=c(36,0),yaxt="n",bty="n",xlim=c(-.2,.7),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)
	
	lines(c(0,0),c(26.5,10.5))
	lines(c(0,0),c(.5,6.5),col="red")
	lines(c(0,0),c(28.5,34.5),col="blue")
	
	axis(side=1,pos=36)
	
	
	
	axis(side=2,at=c(1:7,11:25,29:34),labels=c("SF*w","SB*w","FS*w","FB*w","BS*w","BF*w","w","SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
	                                           "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
	                                           "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
	                                           "SF","SB","FS","FB","BS","BF"),las=1)
	
	temp_95 <- HPDinterval(VCV_with_w_NGM_12t$VCV_Mat[,1:169],prob=.95)
	arrows(temp_95[vect_12t_v2,1],vect_y_v2,temp_95[vect_12t_v2,2],vect_y_v2,code=3,length=.05,angle=90)
	
	temp_80 <- HPDinterval(VCV_with_w_NGM_12t$VCV_Mat[,1:169],prob=.8)
	arrows(temp_80[vect_12t_v2,1],vect_y_v2,temp_80[vect_12t_v2,2],vect_y_v2,code=3,length=0,angle=90,lwd=2,col=vect_col_v2)
	
	points(vectX_12t_v2,vect_y_v2,pch=21,bg="black",cex=.8)
dev.off()


pdf("plots/GZW_12t_only_top.pdf",w=6)
par(mar=c(5,7,4,2))

vect_12t_v2 <- c(157:169)
vectX_12t_v2 <- c(VCV_with_w_NGM_12t$G1_mat)[vect_12t_v2]
vProb <- .95
vect_y_v2 =c(1:6-.2,1:6+.2,7)
vect_col_v2=c(rep("cornflowerblue",6),rep("firebrick3",6),"darkgreen")

plot(vectX_12t_v2,vect_y_v2,ylim=c(7,0),yaxt="n",bty="n",xlim=c(-.2,.7),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)

lines(c(0,0),c(.5,6.5),col="black")


axis(side=1,pos=7.5)

axis(side=2,at=c(1:7),labels=c("SF*w","SB*w","FS*w","FB*w","BS*w","BF*w","w"),las=1)

temp_95 <- HPDinterval(VCV_with_w_NGM_12t$VCV_Mat[,1:169],prob=.95)
arrows(temp_95[vect_12t_v2,1],vect_y_v2,temp_95[vect_12t_v2,2],vect_y_v2,code=3,length=.05,angle=90)

temp_80 <- HPDinterval(VCV_with_w_NGM_12t$VCV_Mat[,1:169],prob=.8)
arrows(temp_80[vect_12t_v2,1],vect_y_v2,temp_80[vect_12t_v2,2],vect_y_v2,code=3,length=0,angle=90,lwd=2,col=vect_col_v2)

points(vectX_12t_v2,vect_y_v2,pch=21,bg="black",cex=.8)
dev.off()

#Now we can compute the Beta estimates - separately for each environment


vect_betas_NGM <- NULL
vect_betas_NaCl <- NULL
k_idx=sample(1:4500,1000)
for(k in 1:1000){
  temp_Gzw <- matrix(VCV_with_w_NGM_12t$VCV_Mat[k_idx[k],1:169],13,13)
  
  vect_betas_NGM <- rbind(vect_betas_NGM ,t(solve(temp_Gzw[1:6,1:6])%*%(temp_Gzw[1:6,13])))	
  vect_betas_NaCl <- rbind(vect_betas_NaCl ,t(solve(temp_Gzw[7:12,7:12])%*%(temp_Gzw[7:12,13])))	  
}

pdf(file='plots/Beta_estimates_all_12_Traits.pdf')
plot(solve(VCV_with_w_NGM_12t$G1_mat[1:12,1:12])%*%(VCV_with_w_NGM_12t$G1_mat[1:12,13]),ylim=c(-1.5,1.5),type="n",ylab=expression(beta),xlab="",xlim=c(0.5,7))

temp_95 <- cbind(apply(vect_betas_NGM,2,function(x){HPDinterval(mcmc(x))}),
                 apply(vect_betas_NaCl,2,function(x){HPDinterval(mcmc(x))}))
temp_80 <- cbind(apply(vect_betas_NGM,2,function(x){HPDinterval(mcmc(x),prob=.8)}),
                 apply(vect_betas_NaCl,2,function(x){HPDinterval(mcmc(x),prob=.8)}))

arrows(c(1:6,1:6+.3),temp_95[1,],c(1:6,1:6+.3), temp_95[2,],code=3,length=.05,angle=90,col='black')
arrows(c(1:6,1:6+.3),temp_80[1,],c(1:6,1:6+.3), temp_80[2,],code=0,length=.0,angle=90,col=rep(c("lightblue","firebrick3"),each=6),lwd=2)

temp_Post <- c(apply(vect_betas_NGM,2,function(x){posterior.mode(mcmc(x))}),
               apply(vect_betas_NaCl,2,function(x){posterior.mode(mcmc(x))}))
points(c(1:6,1:6+.3), temp_Post,pch=16,col='black')

abline(h=0,lty=2)

dev.off()


pdf(file='plots/Beta_estimates_all_12_Traits_only_NaCl.pdf',h=7,w=4)
plot(solve(VCV_with_w_NGM_12t$G1_mat[1:12,1:12])%*%(VCV_with_w_NGM_12t$G1_mat[1:12,13]),ylim=c(-1.5,1.5),type="n",ylab=expression(paste("Selection gradient (", beta[g],")")),xlab="",xlim=c(0.5,6.5),las=1,bty="n",xaxt="n")
axis(side=1,at=1:6,labels=c("SF","SB","FS","FB","BS","BF"))

temp_95 <- cbind(apply(vect_betas_NGM,2,function(x){HPDinterval(mcmc(x))}),
                 apply(vect_betas_NaCl,2,function(x){HPDinterval(mcmc(x))}))
temp_80 <- cbind(apply(vect_betas_NGM,2,function(x){HPDinterval(mcmc(x),prob=.8)}),
                 apply(vect_betas_NaCl,2,function(x){HPDinterval(mcmc(x),prob=.8)}))

arrows(c(1:6),temp_95[1,7:12],c(1:6), temp_95[2,7:12],code=3,length=.05,angle=90,col='black')
arrows(c(1:6),temp_80[1,7:12],c(1:6), temp_80[2,7:12],code=0,length=.0,angle=90,col="firebrick3",lwd=2)

temp_Post <- c(apply(vect_betas_NGM,2,function(x){posterior.mode(mcmc(x))}),
               apply(vect_betas_NaCl,2,function(x){posterior.mode(mcmc(x))}))
points(c(1:6), temp_Post[7:12],pch=16,col='black')

abline(h=0,lty=2)

dev.off()



save(list=ls(),file='output_files/RData/Beta_estimates_all_12t.RData')
