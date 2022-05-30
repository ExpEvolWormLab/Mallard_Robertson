rm(list=ls())
library(MCMCglmm)
library(lme4)
library(data.table)
library(emmeans)
load("output_files/RData/Analysis_Cemee_Pop_WI_NaCl.RData")
v_col = c("chartreuse","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")

########################################################################################
#### We could compute a mean per line estimate in each environment

final_merged_NGM$population2 <- as.character(final_merged_NGM$population)
final_merged_NGM$population2[final_merged_NGM$population2%in%c("GA150","GA250","GA450")] <-  "GA50"
table(final_merged_NGM$population2)

final_merged_NaCl$population2 <- as.character(final_merged_NaCl$population)
final_merged_NaCl$population2[final_merged_NaCl$population2%in%c("CA150","CA250","CA350")] <-  "CA50"
final_merged_NaCl$population2[final_merged_NaCl$population2%in%c("GA150","GA250","GA450")] <-  "GA50"
table(final_merged_NaCl$population2)

final_merged_NGM$population <- as.factor(as.character(final_merged_NGM$population))
final_merged_NGM$population2 <- as.factor(as.character(final_merged_NGM$population2))

final_merged_NaCl$population <- as.factor(as.character(final_merged_NaCl$population))
final_merged_NaCl$population2 <- as.factor(as.character(final_merged_NaCl$population2))


### We can estimate all phenotypic values together / only 2012 effect, specific to the environment ?

sort(names(final_merged))
final_merged=rbind(final_merged_NGM, final_merged_NaCl[,names(final_merged_NaCl)%in%names(final_merged_NGM)])
table(final_merged[,c("location_label","population2")])

final_merged=subset(final_merged,!population2=="CA50")

# Normalize the env. covariates
for(j in c('temperature',"rel_humidity","logD")){
  final_merged[,j] <- (final_merged[,j]-mean(final_merged[,j]))/sd(final_merged[,j])
}

uniqNaCl = unique(subset(final_merged,env_label=='NaCl')$pop_label)
uniqNGM = unique(subset(final_merged,env_label=='NGM')$pop_label)

final_merged_NaCl=subset(final_merged,env_label=='NaCl')
final_merged_NGM=subset(final_merged,env_label=='NGM')

Lines_means=data.frame()
for(i in 1:nb_trait){
print(i)

temp_mod_NGM <- lmer(final_merged_NGM[, vect_P_traits[i]]~ ((temperature+rel_humidity+logD)^3)+is_2012+(1|pop_label)+(1|date_str),data= final_merged_NGM)  
temp_mod_NaCl <- lmer(final_merged_NaCl[, vect_P_traits[i]]~ ((temperature+rel_humidity+logD)^3)+is_2012+(1|pop_label)+(1|date_str),data= final_merged_NaCl)

temp1  <- coef(temp_mod_NGM)[[1]][,1]
temp2  <- row.names(coef(temp_mod_NGM)[[1]])

temp3  <- coef(temp_mod_NaCl)[[1]][,1]
temp4  <- row.names(coef(temp_mod_NaCl)[[1]])

temp_df1=data.frame(pop_label=temp2,trait=as.numeric(temp1))
temp_df2=data.frame(pop_label=temp4,trait=as.numeric(temp3))

names(temp_df1)[2]=paste0(vect_P_traits[i],'_NGM')
names(temp_df2)[2]=paste0(vect_P_traits[i],'_NaCl')

temp_df=merge(temp_df1,temp_df2,all=TRUE)

if(i==1) Lines_means=temp_df
Lines_means <- merge(Lines_means,temp_df,all=TRUE)

}

mean_phen_values_NGM <- Lines_means[,c(1,2,4,6,8,10,12)]
mean_phen_values_NaCl <- Lines_means[,c(1,3,5,7,9,11,13)]


dim(mean_phen_values_NaCl)
dim(mean_phen_values_NGM)

mean_phen_values <- Lines_means
mean_phen_values <- merge(mean_phen_values, unique(final_merged_NaCl[,c("population","pop_label")]))
dim(mean_phen_values)
dim(Lines_means)

mean_phen_values <- mean_phen_values[,c(1,2,4,6,8,10,12,3,5,7,9,11,13,14)]


### Export the phenotypc values and the medians

line_phen_medians <- apply(mean_phen_values[,2:13],2,function(x){
	tapply(x, mean_phen_values$population,median)
})

line_phen_sds <- apply(mean_phen_values[,2:13],2,function(x){
	tapply(x, mean_phen_values$population,sd)
})


save(list=ls(),file='output_files/RData/Phenotypic_analysis_corrected.RData')


##########################################
### Test for ancestral plasticity ?  #####
##########################################

final_merged_NaCl$pop_label <- paste(final_merged_NaCl$pop_label,"NaCl",sep="_")
final_merged_NGM$pop_label <- paste(final_merged_NGM$pop_label,"NGM",sep="_")

selected_col <- c(vect_P_traits,"population", "population2","temperature","rel_humidity","logD","pop_label","date_str","env_label","is_2012")
final_export_all <- rbind(final_merged_NaCl[,selected_col],final_merged_NGM[,selected_col])

head(final_export_all)

final_export_all$temperature= (final_export_all$temperature-mean(final_export_all $temperature))/sd(final_export_all $temperature)
final_export_all $rel_humidity= (final_export_all $rel_humidity-mean(final_export_all $rel_humidity))/sd(final_export_all $rel_humidity)
final_export_all $logD= (final_export_all $logD-mean(final_export_all$logD))/sd(final_export_all $logD)

final_export_all_A6140 <- subset(final_export_all,population=="A6140")
dim(final_export_all_A6140)
delta_plasticity_lmer <- NULL
plast_estimates_for_plot =NULL
for(i in 1:6){
temp_mod1 <- lmer(final_export_all_A6140[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + env_label + is_2012 + (1|pop_label) +(1|date_str),data= final_export_all_A6140)
print(summary(temp_mod1)$coef[1,1])
temp_mod2 <- lmer(final_export_all_A6140[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + is_2012 + (1|pop_label) +(1|date_str),data= final_export_all_A6140)
print(anova(temp_mod1,temp_mod2))

delta_plasticity_lmer <- c(delta_plasticity_lmer,-as.numeric(coef(temp_mod1)[[1]][1,][names(coef(temp_mod1)[[1]][1,])=="env_labelNGM"]))

# 95% CI for the plot
temp_alt= lmer(final_export_all_A6140[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + env_label + is_2012 + (1|pop_label) + (1|date_str) -1 ,data= final_export_all_A6140)
plast_estimates_for_plot=rbind(plast_estimates_for_plot,cbind(summary(temp_alt)$coef[4:5,1],confint(temp_alt)[7:8,]))

}

final_export_all$pop_label_init <- tstrsplit(final_export_all$pop_label,"_")[[1]]
final_export_all_NaCl <- subset(final_export_all,env_label=='NaCl')

temp_mod1 <- lmer(T12~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NaCl)
emT12 <- emmeans(temp_mod1,~population)
pairs(emT12)

temp_mod1 <- lmer(T13~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NaCl)
emT13 <- emmeans(temp_mod1,~population)
pairs(emT13)


temp_mod1 <- lmer(T21~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NaCl)
emT21 <- emmeans(temp_mod1,~population)
pairs(emT21)


temp_mod1 <- lmer(T23~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NaCl)
emT23 <- emmeans(temp_mod1,~population)
pairs(emT23)


temp_mod1 <- lmer(T31~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NaCl)
emT31 <- emmeans(temp_mod1,~population)
pairs(emT31)


temp_mod1 <- lmer(T32~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NaCl)
emT32 <- emmeans(temp_mod1,~population)
pairs(emT32)


### Final Figures

pdf("plots/Phenotypic_plasticity.pdf")

par(mfrow=c(2,3))
for(k in 1:6){
  plot(rep(0,nrow(mean_phen_values)),mean_phen_values[,k+1],xlim=c(-1,2),type="n",bty="n",las=1,xaxt="n",xlab="",ylab="log transition rates")
  axis(side=1,at=c(0,1),labels=c("Low\nSalt","High\nSalt"),padj=.5)
  temp <- subset(mean_phen_values,population=="A6140")
  points(jitter(rep(0,nrow(temp))),temp[,k+1],pch=16,cex=.8,col="gray")
  points(1+jitter(rep(0,nrow(temp))),temp[,k+1+6],pch=16,cex=.8,col="gray")
  
  arrows(0 ,plast_estimates_for_plot[(2*k),2], 0 ,plast_estimates_for_plot[(2*k),3],code=3,angle=90,length=.1,col='black')
  arrows(1 ,plast_estimates_for_plot[(2*k-1),2], 1 ,plast_estimates_for_plot[(2*k-1),3],code=3,angle=90,length=.1,col='black')
  points(c(0,1),plast_estimates_for_plot[(2*k):(2*k-1),1],pch=21,bg="orange",cex=1,col="black",type="b",lwd=2)
  
}

dev.off()


### And now the estimation of trait change during evolution in NaCl

delta_evolution_lmer <- NULL
delta_evolution_lmer_GA <- NULL

evol_estimates_for_plot=NULL
evol_estimates_for_plot_GA_all=NULL

for(i in 1:6){

temp_mod1 <- lmer(final_export_all_NaCl[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NaCl)
temp_mod2 <- lmer(final_export_all_NaCl[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + is_2012 + (1|pop_label) +(1|date_str),data= final_export_all_NaCl)
print(anova(temp_mod1,temp_mod2))

delta_evolution_lmer_GA <- rbind(delta_evolution_lmer_GA,as.numeric(coef(temp_mod1)[[1]][1,][names(coef(temp_mod1)[[1]][1,])%in%c("populationGA150","populationGA250","populationGA450")]))
# 95% CI for the plot
temp_alt= lmer(final_export_all_NaCl[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + population + is_2012   + (1|pop_label) +(1|date_str)-1,data= final_export_all_NaCl)
evol_estimates_for_plot_GA_all=rbind(evol_estimates_for_plot_GA_all,cbind((data.frame(summary(temp_alt)$coef[4:7,1])),confint(temp_alt)[7:10,]))

temp_mod1 <- lmer(final_export_all_NaCl[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + is_2012 + population2 + (1|pop_label) +(1|date_str)+(1|population),data= final_export_all_NaCl)
temp_mod2 <- lmer(final_export_all_NaCl[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + is_2012 + 1 + (1|pop_label) +(1|date_str)+(1|population),data= final_export_all_NaCl)
print(anova(temp_mod1,temp_mod2))

# 95% CI for the plot
temp_alt= lmer(final_export_all_NaCl[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + population2 -1 + is_2012   + (1|pop_label) +(1|date_str)+(1|population),data= final_export_all_NaCl)
evol_estimates_for_plot=rbind(evol_estimates_for_plot,cbind(summary(temp_alt)$coef[4:5,1],confint(temp_alt)[8:9,]))

delta_evolution_lmer <- c(delta_evolution_lmer,as.numeric(coef(temp_mod1)[[1]][1,][names(coef(temp_mod1)[[1]][1,])%in%c("population2GA50")]))

}


############ Phenotypic Evolution

pdf("plots/Phenotypic_evolution_NaCl.pdf")

par(mfrow=c(2,3))
for(k in 1:6){

  plot(rep(0,nrow(mean_phen_values)),mean_phen_values[,k+1],xlim=c(-1,2),type="n",bty="n",las=1,xaxt="n",xlab="",ylab="")
  axis(side=1,at=c(0,1),labels=c("A6140","GA"),padj=.5)
  temp <- subset(mean_phen_values,population=="A6140")

  points(rep(0,nrow(temp)),temp[,k+1+6],pch=16,cex=.8,col="gray")
  
  for(p in 1:3){
    temp <- subset(mean_phen_values,population==c("GA150","GA250","GA450")[p])	
    points(rep(1+.1+0.1*p,nrow(temp)),temp[,k+1+6],pch=16,cex=.8,col=v_col[p+4])
  }
  
  
  arrows(.9 ,evol_estimates_for_plot[(2*k),2], .9 ,evol_estimates_for_plot[(2*k),3],code=3,angle=90,length=.08,col='black')
  arrows(0 ,evol_estimates_for_plot[(2*k-1),2], 0 ,evol_estimates_for_plot[(2*k-1),3],code=3,angle=90,length=.08,col='black')
  
  points(c(1.1+0.1*(1:3)),evol_estimates_for_plot_GA_all[(4*k-2):(4*k),1],bg=v_col[5:7],pch=21)
  points(c(0,.9),evol_estimates_for_plot[(2*k-1):(2*k),1],pch=16,type="b")
  

}

dev.off()


### And now the estimation of trait change during evolution in NGM

final_export_all_NGM <- subset(final_export_all,env_label=='NGM')

delta_evolution_lmer_NGM <- NULL
delta_evolution_lmer_GA_NGM <- NULL

evol_estimates_for_plot_NGM=NULL
evol_estimates_for_plot_GA_all_NGM=NULL

for(i in 1:6){
  
  temp_mod1 <- lmer(final_export_all_NGM[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NGM)
  temp_mod2 <- lmer(final_export_all_NGM[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + is_2012 + (1|pop_label) +(1|date_str),data= final_export_all_NGM)
  print(anova(temp_mod1,temp_mod2))
  
  # 95% CI for the plot
  temp_alt= lmer(final_export_all_NGM[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + population + is_2012   + (1|pop_label) +(1|date_str)-1,data= final_export_all_NGM)
  evol_estimates_for_plot_GA_all_NGM=rbind(evol_estimates_for_plot_GA_all_NGM,cbind((data.frame(summary(temp_alt)$coef[4:7,1])),confint(temp_alt)[7:10,]))
  
  delta_evolution_lmer_GA_NGM <- rbind(delta_evolution_lmer_GA_NGM,as.numeric(coef(temp_mod1)[[1]][1,][names(coef(temp_mod1)[[1]][1,])%in%c("populationGA150","populationGA250","populationGA450")]))
  
  temp_mod1 <- lmer(final_export_all_NGM[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + is_2012 + population2 + (1|pop_label) +(1|date_str)+(1|population),data= final_export_all_NGM)
  temp_mod2 <- lmer(final_export_all_NGM[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + is_2012 + 1 + (1|pop_label) +(1|date_str)+(1|population),data= final_export_all_NGM)
  print(anova(temp_mod1,temp_mod2))
  
  # 95% CI for the plot
  temp_alt= lmer(final_export_all_NGM[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + population2 -1 + is_2012   + (1|pop_label) +(1|date_str)+(1|population),data= final_export_all_NGM)
  evol_estimates_for_plot_NGM=rbind(evol_estimates_for_plot_NGM,cbind(summary(temp_alt)$coef[4:5,1],confint(temp_alt)[8:9,]))
  
  
  delta_evolution_lmer_NGM <- c(delta_evolution_lmer_NGM,as.numeric(coef(temp_mod1)[[1]][1,][names(coef(temp_mod1)[[1]][1,])%in%c("population2GA50")]))
  
  
  
}



pdf("plots/Phenotypic_evolution_NGM.pdf")

par(mfrow=c(2,3))
for(k in 1:6){
  
  plot(rep(0,nrow(mean_phen_values)),mean_phen_values[,k+1],xlim=c(-1,2),type="n",bty="n",las=1,xaxt="n",xlab="",ylab="")
  axis(side=1,at=c(0,1),labels=c("A6140","GA"),padj=.5)
  temp <- subset(mean_phen_values,population=="A6140")
  
  points(rep(0,nrow(temp)),temp[,k+1],pch=16,cex=.8,col="gray")
  
  
  for(p in 1:3){
    temp <- subset(mean_phen_values,population==c("GA150","GA250","GA450")[p])	
    points(rep(1+.1+0.1*p,nrow(temp)),temp[,k+1],pch=16,cex=.8,col=v_col[p+4])
  }
  
  
  arrows(.9 ,evol_estimates_for_plot_NGM[(2*k),2], .9 ,evol_estimates_for_plot_NGM[(2*k),3],code=3,angle=90,length=.08,col='black')
  arrows(0 ,evol_estimates_for_plot_NGM[(2*k-1),2], 0 ,evol_estimates_for_plot_NGM[(2*k-1),3],code=3,angle=90,length=.08,col='black')
  
  points(c(1.1+0.1*(1:3)),evol_estimates_for_plot_GA_all_NGM[(4*k-2):(4*k),1],bg=v_col[5:7],pch=21)
  points(c(0,.9),evol_estimates_for_plot_NGM[(2*k-1):(2*k),1],pch=16,type="b")
  
  
}

dev.off()

temp_mod1 <- lmer(T12~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NGM)
emT12 <- emmeans(temp_mod1,~population)
pairs(emT12)

temp_mod1 <- lmer(T13~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NGM)
emT13 <- emmeans(temp_mod1,~population)
pairs(emT13)


temp_mod1 <- lmer(T21~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NGM)
emT21 <- emmeans(temp_mod1,~population)
pairs(emT21)


temp_mod1 <- lmer(T23~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NGM)
emT23 <- emmeans(temp_mod1,~population)
pairs(emT23)


temp_mod1 <- lmer(T31~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NGM)
emT31 <- emmeans(temp_mod1,~population)
pairs(emT31)


temp_mod1 <- lmer(T32~ (temperature+rel_humidity+logD)^3 + is_2012 + population + (1|pop_label) +(1|date_str),data= final_export_all_NGM)
emT32 <- emmeans(temp_mod1,~population)
pairs(emT32)



save(list=c("final_export_all","line_phen_medians","line_phen_sds","mean_phen_values","delta_plasticity_lmer","delta_evolution_lmer_GA","delta_evolution_lmer","delta_evolution_lmer_GA_NGM","delta_evolution_lmer_NGM"),file="output_files/RData/Analysis_Cemee_Pop_WI_Line_Means_corrected.RData")
