rm(list=ls());gc()
library(matrixStats)

load("output_files/RData/Phenotypic_data.RData")
vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")

env_covariates = c("temperature","logD","rel_humidity","env_label")
retained_names=c(vect_P_traits,env_covariates,'date_str',"pop_label","is_2012","population")

all_data <- final_merged[,names(final_merged)%in%retained_names]
dim(all_data) # 1366
all_data=na.omit(all_data)
dim(all_data) # 1366

all_data$line_env = paste(all_data$pop_label,all_data$env_label,sep='_')
all_data$line_pop = paste(all_data$pop_label,all_data$population,sep='_')
all_data$pop_env = paste(all_data$population,all_data$env_label,sep='_')

all_data$date_str = as.factor(as.character(all_data$date_str))
all_data$env_label = factor(all_data$env_label,levels=c("NaCl","NGM"))
##########

nb_trait = length(vect_P_traits)

for(j in c('temperature',"rel_humidity","logD")){
  print(j)
  all_data[,j] <- (all_data[,j]-mean(all_data[,j]))/sd(all_data[,j])
}

all_data_A6140=(subset(all_data,population=="A6140"))
unique.manova <- manova(cbind(T12, T13, T21, T23, T31, T32,area.F) ~ env_label*population + (temperature+rel_humidity+logD)^3 + is_2012 + date_str + Error(line_env) , data=all_data)
summary(unique.manova)

## Save the coefficients in a table
write.table(unique.manova$line_env$coefficients,file='output_files/txt/Manova.estimates.txt',sep='\t',quote=FALSE)
write.table(capture.output(summary.aov(unique.manova$line_env)),file='output_files/txt/Manova_trait_responses.txt',row.names=FALSE,quote=FALSE,sep='\t')

# summary for the block effects
id_block = substring(row.names(unique.manova$line_env$coefficients),1,4)=="date"
colMeans(unique.manova$line_env$coefficients[id_block,],na.rm=TRUE)
colSds(unique.manova$line_env$coefficients[id_block,],na.rm=TRUE)

m.unique = summary(unique.manova, test = "Wilks")
write.csv(capture.output(m.unique$`Error: line_env`),file='output_files/txt/Manova_results.txt',quote=FALSE,row.names=FALSE)

# Calculate Dplast as the SSCP matrix for env. effect
Dplast <- m.unique$`Error: line_env`$SS$env_label
Ddiv <- m.unique$`Error: line_env`$SS$population

EV_plast <- eigen(Dplast)$vectors
EV_div <- eigen(Ddiv)$vectors # oriented to get an increase of both divergence and plastic traits

all_data$dplast <- t(EV_plast[,1] %*% t(all_data[,vect_P_traits]))
all_data$ddiv <- t(EV_div[,1] %*% t(all_data[,vect_P_traits]))

(eigen(Dplast)$values[1])/sum(eigen(Dplast)$values) # 1
eigen(Ddiv)$values[1]/sum(eigen(Ddiv)$values); eigen(Ddiv)$values[2]/sum(eigen(Ddiv)$values)

save(list=c("EV_div","EV_plast"),file='output_files/RData/SSCP_eigenvectors.RData')

write.table(file='output_files/txt/SSCP_plasticity.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE,round(Dplast,digits=3))
write.table(file='output_files/txt/SSCP_divergence.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE,round(Ddiv,digits=3))

write.table(file='output_files/txt/SSCP_plasticity_ED.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE,round(EV_plast,digits=3))
write.table(file='output_files/txt/SSCP_divergence_ED.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE,round(EV_div,digits=3))

save(list=c("all_data"),file='output_files/RData/Phenotypic_data_for_MANOVA.RData')
