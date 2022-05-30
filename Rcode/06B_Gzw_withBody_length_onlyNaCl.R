rm(list=ls());gc()
library(MCMCglmm)
library(lme4)
library(data.table)
library(matrixStats)


## Additinal code to look at Gzw matrix with body size - not shown in the manuscript

load('output_files/RData/Analysis_Cemee_Pop_WI_NaCl.RData')
merged_Phen=read.table("data/merged_data_means.txt",h=TRUE)

shared_names=names(final_merged)[names(final_merged)%in%names(merged_Phen)]

temp.merged <- (merge(final_merged,merged_Phen,all.x=TRUE))
temp.merged.missing <- subset(temp.merged,is.na(curvature.F.var))[,c("exper_data_id",shared_names)]

data_for_G_NaCl <- (merge(data_for_G_NaCl,merged_Phen,all.x=TRUE))
data_for_G_NGM <- (merge(data_for_G_NGM,merged_Phen,all.x=TRUE))

load("data/fertility.rda")

mean_relative_fit <- subset(ecoefs, pop=="A6140" & env=="NaCl")
mean_relative_fit$w <- exp(mean_relative_fit$fertility)/mean(exp(mean_relative_fit$fertility))

data_for_G_NGM_A6140 <- subset(data_for_G_NGM,population=="A6140")
data_for_G_NaCl_A6140 <- subset(data_for_G_NaCl,population=="A6140")

names(data_for_G_NGM_A6140)[names(data_for_G_NGM_A6140)%in%vect_P_traits] <- paste0(vect_P_traits,"_NGM")
names(data_for_G_NaCl_A6140)[names(data_for_G_NaCl_A6140)%in%vect_P_traits] <- paste0(vect_P_traits,"_NaCl")

names(data_for_G_NGM_A6140)[names(data_for_G_NGM_A6140)=='area.F']='area.F_NGM'
names(data_for_G_NaCl_A6140)[names(data_for_G_NaCl_A6140)=='area.F']='area.F_NaCl'

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
vect_P_traits_all = paste0(rep(vect_P_traits,1),rep(c("_NaCl"),each=6)) 

### area.F

retained_names=c(vect_P_traits_all,env_covariates,'date_str',"pop_label","is_2012","area.F_NGM","area.F_NaCl")

all_data <- data_for_G_NaCl_A6140[,names(data_for_G_NaCl_A6140)%in%retained_names]
dim(all_data) # 357 13

mean_relative_fit <- mean_relative_fit[,c("line","w")]
names(mean_relative_fit)[1] = 'pop_label'

dim(all_data)
all_data_with_w <- (merge(all_data, mean_relative_fit))
dim(all_data_with_w) # 334 14

all_data_with_w=na.omit(all_data_with_w)
dim(all_data_with_w) # 328 24

vect_P_traits_all=c(vect_P_traits_all,"area.F_NaCl","w")
length(vect_P_traits_all)
nb_trait = length(vect_P_traits_all)

### Normalize ?

meanSd_TR <- mean(colSds(as.matrix(all_data_with_w[,vect_P_traits_all[1:6]])))
sd(as.matrix(all_data_with_w[,vect_P_traits_all[7]]))

for(i in 1:6){
  all_data_with_w[,vect_P_traits_all[i]] <- (all_data_with_w[,vect_P_traits_all[i]]-mean(all_data_with_w[,vect_P_traits_all[i]]))/meanSd_TR
}

mean(colSds(as.matrix(all_data_with_w[,vect_P_traits_all[1:6]])))

i=7
all_data_with_w[,vect_P_traits_all[i]] <- (all_data_with_w[,vect_P_traits_all[i]]-mean(all_data_with_w[,vect_P_traits_all[i]],na.rm=TRUE))/sd(as.matrix(all_data_with_w[,vect_P_traits_all[7]]))
colSds(as.matrix(all_data_with_w[,vect_P_traits_all]))

k=0

dim(all_data_with_w[,vect_P_traits_all])
	phen.var = diag(nb_trait) * diag(var(subset(all_data_with_w, select = vect_P_traits_all)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC <- MCMCglmm(cbind(c(T12_NaCl, T13_NaCl, T21_NaCl, T23_NaCl, T31_NaCl, T32_NaCl,area.F_NaCl,w)) ~ 
	                          (temperature_NaCl+rel_humidity_NaCl+logD_NaCl)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = all_data_with_w, prior = prior_mod, verbose = TRUE,nitt=1000000, burnin=100000,thin=200)

	save(list=ls(),file='output_files/RData/Beta_NaCl_withBL_scaled.RData')
	
	post_dist = posterior.mode(model_MCMC$VCV)

	VCV_with_w_NaCl_with_BL <- list(Population = "A6_with_12t_with_BL", N_measurement = nrow(temp_final), G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)

	VCV_with_w_NaCl_with_BL$G1_mat


pdf(file='plots/GZW_NaCl_with_BL_scaled.pdf',h=9,w=5.5)
	
	par(mar=c(5,7,4,2))
	
	vect_6t_v2_with_BL <- c(c(2:7,11:15,20:23,29:31,38,39,47),c(57:63),c(9*(1:8)-8))
	vectX_6t_v2 <- c(VCV_with_w_NaCl_with_BL$G1_mat)[vect_6t_v2_with_BL]
	
	vProb <- .95
	
	vect_y_v2 =c(c(11:31),c(1:7),c(35:41),8)
	vect_col_v2=(c(rep("lightblue",21),rep("darkgreen",7),rep("lightblue",7),"darkgreen"))
	
	length(vect_y_v2)
	length(vect_col_v2)

	plot(vect_6t_v2_with_BL,vect_y_v2,ylim=c(44,0),yaxt="n",bty="n",xlim=c(-.4,.7),xlab="Genetic (co-)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)
	mtext(side=2,"Phenotypic traits\n Diagonal                                    Off-diagonal                                       Fitness      ",padj=-2,cex=1.2)
	lines(c(0,0),c(31.5,10.5))
	lines(c(0,0),c(.5,7.5),col="red")
	lines(c(0,0),c(34.5,41.5),col="blue")
	
	axis(side=1,pos=44)
	
	axis(side=2,at=c(1:8,11:31,35:41),labels=c("SF*w","SB*w","FS*w","FB*w","BS*w","BF*w","Area*w","w",
	                                           "SF*SB","SF*FS","SF*FB","SF*BS","SF*BF","SF*Area",
	                                           "SB*FS","SB*FB","SB*BS","SB*BF","SB*Area",
	                                           "FS*FB","FS*BS","FS*BF","FS*Area",
	                                           "FB*BS","FB*BF","FB*Area","BS*BF","BS*Area","BF*Area",
	                                           "SF","SB","FS","FB","BS","BF","Area"),las=1)
	
	temp_95 <- HPDinterval(VCV_with_w_NaCl_with_BL$VCV_Mat[,1:nb_trait^2],prob=.95)
	temp_80 <- HPDinterval(VCV_with_w_NaCl_with_BL$VCV_Mat[,1:nb_trait^2],prob=.8)
	
	arrows(temp_95[vect_6t_v2_with_BL,1],vect_y_v2,temp_95[vect_6t_v2_with_BL,2],vect_y_v2,code=3,length=.05,angle=90)
	arrows(temp_80[vect_6t_v2_with_BL,1],vect_y_v2,temp_80[vect_6t_v2_with_BL,2],vect_y_v2,code=3,length=0,angle=90,lwd=2,col=vect_col_v2)
	points(vectX_6t_v2,vect_y_v2,pch=21,bg="black",cex=.8)
	
	
	
	dev.off()
	
	save(list=ls(),file='output_files/RData/Beta_NaCl_withBL_scaled.RData')
	
	#Now we can compute the Beta estimates - separately for each environment
	
	vect_betas_withBL <- NULL
	k_idx=sample(1:4500,1000)
	for(k in 1:1000){
	  temp_Gzw <- matrix(VCV_with_w_NaCl_with_BL$VCV_Mat[k_idx[k],1:64],8,8)
	  vect_betas_withBL <- rbind(vect_betas_withBL ,t(solve(temp_Gzw[1:7,1:7])%*%(temp_Gzw[1:7 ,8])))	
	}
	
	### Or if we compute the beta estimates with the former G ?
	
	pdf(file='plots/Beta_estimates_NaCl_withBL.pdf')
	plot(solve(VCV_with_w_NaCl_with_BL$G1_mat[1:7,1:7])%*%(VCV_with_w_NaCl_with_BL$G1_mat[1:7,8]),ylim=c(-1.5,1.5),type="n",ylab=expression(beta),xlab="",xlim=c(0.5,7))
	
	temp_95 <- apply(vect_betas_withBL,2,function(x){HPDinterval(mcmc(x))})
	temp_80 <- apply(vect_betas_withBL,2,function(x){HPDinterval(mcmc(x),prob=.8)})
	            
	arrows(c(1:7),temp_95[1,],c(1:7), temp_95[2,],code=3,length=.05,angle=90,col='black')
	arrows(c(1:7),temp_80[1,],c(1:7), temp_80[2,],code=0,length=.0,angle=90,col=c(rep("lightblue",6),"yellow"),lwd=2)
	
	temp_Post <- apply(vect_betas_withBL,2,function(x){posterior.mode(mcmc(x))})
	points(1:7, temp_Post,pch=16,col='black')
	
	abline(h=0,lty=2)
	
	dev.off()
	
	