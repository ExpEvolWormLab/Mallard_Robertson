rm(list=ls());gc()
library(MCMCglmm)
library(lme4)
library(data.table)
library(emmeans)

## Additinal code to look at body size evolution - not shown in the manuscript

load('output_files/RData/Analysis_Cemee_Pop_WI_NaCl.RData')
# Load all other phenotypes
merged_Phen=read.table("data/merged_data_means.txt",h=TRUE)

merged_Phen=subset(merged_Phen,substring(pop_label,1,2)%in%c("A6","GA") & env_label=='NaCl' & data_group_name!='B_maleCtrl')
dim(merged_Phen) # 728
merged_Phen$population2=substring(merged_Phen$pop_label,1,2)
merged_Phen$population=substring(merged_Phen$pop_label,1,5)

table(substring(merged_Phen$date_str,1,4))
merged_Phen$is_2012= (substring(merged_Phen$date_str,1,4)=="2012")
merged_Phen=subset(merged_Phen,substring(merged_Phen$date_str,1,4)!="2016")

for(j in c('temperature',"rel_humidity","logMeanTracksPerSecond")){
  merged_Phen[,j] <- (merged_Phen[,j]-mean(merged_Phen[,j]))/sd(merged_Phen[,j])
}

temp_mod1 <- lmer(area.F~ (logMeanTracksPerSecond + temperature + rel_humidity) + is_2012  + population2 + (1|pop_label) +(1|date_str)+(1|population),data= merged_Phen)
temp_mod2 <- lmer(area.F~ (logMeanTracksPerSecond + temperature + rel_humidity) + is_2012  + (1|pop_label) +(1|date_str)+(1|population),data= merged_Phen)
anova(temp_mod1,temp_mod2)
drop1(temp_mod1)

##########

temp_mod1 <- lmer(area.F~ (logMeanTracksPerSecond + temperature + rel_humidity)^3 + is_2012  + population + (1|pop_label) +(1|date_str),data= merged_Phen)
temp_mod2 <- lmer(area.F~ (logMeanTracksPerSecond + temperature + rel_humidity)^3 + is_2012  + (1|pop_label) +(1|date_str),data= merged_Phen)
anova(temp_mod1,temp_mod2)
drop1(temp_mod1)

summary(temp_mod1)

em_area <- emmeans(temp_mod1,~population)
pairs(em_area)
pdf("plots/Raw_Area_per_population.pdf")
boxplot(merged_Phen$area.F~ merged_Phen$population)
dev.off()

#### One is significantly different
temp_mod1 <- lmer(area.F~ (logMeanTracksPerSecond + temperature + rel_humidity)^3 + is_2012  + population2 +(1|population) + (1|pop_label) +(1|date_str),data= merged_Phen)
summary(temp_mod1)

## FOr the manuscript : 
#Standardized phen. evolution
(2.447e-03)/sd(merged_Phen$area.F)

