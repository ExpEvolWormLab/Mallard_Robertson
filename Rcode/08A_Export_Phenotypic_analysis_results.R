rm(list=ls());gc()
library(lme4)
library(emmeans)

## This code produce plots for figures 1 and 6 and export some reusults tables 

load('output_files/RData/Phenotypic_data_for_MANOVA.RData')
vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")

model_line_means = unique(all_data[,c("pop_label","env_label","line_env")])
model_line_means=merge(model_line_means,unique(all_data[,c("pop_label","population")]))
dim(model_line_means) # 703

model_line_means$temperature=0
model_line_means$rel_humidity=0
model_line_means$logD=0
model_line_means$is_2012=TRUE

model_populations_means = unique(model_line_means[,c("env_label","population")])
model_populations_plasticity = data.frame(population=unique(model_line_means[,'population']))
plast_constrast=NULL
divergence_constrast=NULL
divergence_constrast_NGM=NULL

# Order the population means to add the estimates with CI
model_populations_means = model_populations_means[order(model_populations_means$env_label,model_populations_means$population),]

for(i in 1:7){
trait_data <- all_data[,vect_P_traits[i]]
mod_temp = lmer( trait_data ~ population*env_label + (temperature+rel_humidity+logD)^3 + is_2012 + (1|date_str) + (1|line_env),data= all_data)
#Store appropriate comparisons
em_population <- as.data.frame(pairs(emmeans(mod_temp,~population*env_label)))
plast_constrast = rbind(plast_constrast,subset(em_population,contrast=='A6140 NaCl - A6140 NGM'))
divergence_constrast = rbind(divergence_constrast,subset(em_population,contrast%in%paste0('A6140 NaCl - GA',c(1,2,4),'50 NaCl')))
divergence_constrast_NGM = rbind(divergence_constrast_NGM,subset(em_population,contrast%in%paste0('A6140 NGM - GA',c(1,2,4),'50 NGM')))

#Store the Least-square populations estimates
temp_line_means = as.data.frame(emmeans(mod_temp,specs=c("population","env_label")))
names(temp_line_means)[c(3,6,7)] = paste0(vect_P_traits[i],"_",names(temp_line_means)[c(3,6,7)])
model_populations_means=cbind(model_populations_means,temp_line_means[,c(3,6,7)])

temp_line_means_plasticity = as.data.frame(emmeans(mod_temp,specs=c("population","env_label")))
emmeans(mod_temp,specs=c("population","env_label"))

confint(contrast(emmeans(mod_temp,~population:env_label),'pairwise',by='population'))
# Get the line_env BLUPs
model_line_means=merge(model_line_means,
data.frame(line_env=row.names(coef(mod_temp)[[1]]),temp=coef(mod_temp)[[1]][,1]))
names(model_line_means)[ncol(model_line_means)]=vect_P_traits[i]

}

# Add the emmeans effects to the BLUPs (and remove the intercept from the lmer model)
for(i in 1:nrow(model_populations_means)){
  is_select <- model_line_means$env_label==model_populations_means$env_label[i] & model_line_means$population==model_populations_means$population[i]
 for(k in 1:7){
   model_line_means[is_select,vect_P_traits[k]] = model_line_means[is_select,vect_P_traits[k]] - mean(model_line_means[is_select,vect_P_traits[k]]) + model_populations_means[i,paste0(vect_P_traits[k],"_emmean")]
 }
}

for(i in c(2:5)) plast_constrast[,i]=round(plast_constrast[,i],digits=3)
for(i in c(2:5)) divergence_constrast[,i]=round(divergence_constrast[,i],digits=3)
divergence_constrast_NGM[,2:5]=round(divergence_constrast_NGM[,2:5],digits=3)

## Export contrast tables:
write.table(plast_constrast,file='output_files/txt/Plasticity_contrasts.txt',sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(divergence_constrast,file='output_files/txt/Divergence_contrasts_High_Salt.txt',sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(divergence_constrast_NGM,file='output_files/txt/Divergence_contrasts_Low_Salt.txt',sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

# Plot for Figure 1
pdf("plots/Phenotypic_plasticity.pdf",h=10)
par(mfrow=c(3,3),mar=c(5,5,4,2))
for(k in 1:7){
  
  if(k!=7){
    plot(rep(0,nrow(model_line_means)),model_line_means[,vect_P_traits[k]],xlim=c(-1,2),type="n",bty="n",las=1,xaxt="n",xlab="",ylab="Transition rates (log Hz)",
         main=c("SF","SB","FS","FB","BS",'BF','Size')[k])
  }else{
    plot(rep(0,nrow(model_line_means)),model_line_means[,vect_P_traits[k]],xlim=c(-1,2),type="n",bty="n",las=1,xaxt="n",xlab="",ylab=expression(paste("Area  (50 x ",mm^2,")")),
         main=c("SF","SB","FS","FB","BS",'BF','Size')[k])
  }
  
  #if(k>3) axis(side=1,at=c(0,1),labels=c("Low\nSalt","High\nSalt"),padj=.5)
  axis(side=1,at=c(0,1),labels=c("Low\nSalt","High\nSalt"),padj=.5)
  
  temp <- subset(model_line_means,population=="A6140" & env_label=="NGM")
  temp2 <- subset(model_line_means,population=="A6140" & env_label=="NaCl")

  points(jitter(c(rep(0,nrow(temp)),rep(1,nrow(temp2)))),rbind(temp,temp2)[,vect_P_traits[k]],pch=16,cex=.8,col="gray")
  
  arrows(c(1,0) , model_populations_means[model_populations_means$population=="A6140",paste0(vect_P_traits[k],"_lower.CL")], c(1,0) ,model_populations_means[model_populations_means$population=="A6140",paste0(vect_P_traits[k],"_upper.CL")],code=3,angle=90,length=.1,col='black')
  if(k!=4){
    points(c(1,0),model_populations_means[model_populations_means$population=="A6140",paste0(vect_P_traits[k],"_emmean")],pch=21,bg="grey",cex=1,col="black",type="b",lwd=2)
  }else{
    points(c(1,0),model_populations_means[model_populations_means$population=="A6140",paste0(vect_P_traits[k],"_emmean")],pch=21,bg="grey",cex=1,col="black",type="p",lwd=2)
  }
  
}
dev.off()


col_GA <- c("cadetblue1", "cornflowerblue", "slateblue2")

#### Now divergence

# Plot for Figure 6
pdf("plots/Phenotypic_divergence_High_Salt.pdf",w=8)
par(mfrow=c(3,3))

posX = c(0,.7,1,1.3)

for(k in 1:7){
  
  model_line_means_NaCl <- subset(model_line_means, env_label=='NaCl')
  
  if(k!=7){
    plot(rep(0,nrow(model_line_means_NaCl)),model_line_means_NaCl[,vect_P_traits[k]],xlim=c(-1,2),type="n",bty="n",las=1,xaxt="n",xlab="",ylab="Transition rates (log Hz)",
         main=c("SF","SB","FS","FB","BS",'BF','Size')[k])
  }else{
    plot(rep(0,nrow(model_line_means_NaCl)),model_line_means_NaCl[,vect_P_traits[k]],xlim=c(-1,2),type="n",bty="n",las=1,xaxt="n",xlab="",ylab=expression(paste("Area  (50 x ",mm^2,")")),
         main=c("SF","SB","FS","FB","BS",'BF','Size')[k])
  }
  

  axis(side=1,at=c(0,1),labels=c("A6140","GA[1,2,4]50"),padj=.5)
  
  temp <- subset(model_line_means,population=="A6140" & env_label=='NaCl')
  temp2 <- subset(model_line_means,population=="GA150" & env_label=='NaCl')
  temp3 <- subset(model_line_means,population=="GA250" & env_label=='NaCl')
  temp4 <- subset(model_line_means,population=="GA450" & env_label=='NaCl')
  
  points(jitter(rep(posX,c(nrow(temp),nrow(temp2),nrow(temp3),nrow(temp4)))),rbind(temp,temp2,temp3,temp4)[,vect_P_traits[k]],
         pch=16,cex=.8,col=rep(c("grey",col_GA),c(nrow(temp),nrow(temp2),nrow(temp3),nrow(temp4))))
  
  
  arrows(posX , model_populations_means[model_populations_means$env_label=="NaCl",paste0(vect_P_traits[k],"_lower.CL")], posX,model_populations_means[model_populations_means$env_label=="NaCl",paste0(vect_P_traits[k],"_upper.CL")],code=3,angle=90,length=.1,col='black')
  points(posX , model_populations_means[model_populations_means$env_label=="NaCl",paste0(vect_P_traits[k],"_emmean")],pch=21,bg="white",cex=1,col="black",type="b",lwd=2)
  
}
dev.off()


# Plot for Figure 6 - figure supplement 2
pdf("plots/Phenotypic_divergence_Low_Salt.pdf",w=8)
par(mfrow=c(3,3))
posX = c(0,.7,1,1.3)

for(k in 1:7){
  model_line_means_NGM <- subset(model_line_means, env_label=='NGM')
  
  if(k!=7){
    plot(rep(0,nrow(model_line_means_NGM)),model_line_means_NGM[,vect_P_traits[k]],xlim=c(-1,2),type="n",bty="n",las=1,xaxt="n",xlab="",ylab="Transition rates (log Hz)",
         main=c("SF","SB","FS","FB","BS",'BF','Size')[k])
  }else{
    plot(rep(0,nrow(model_line_means_NGM)),model_line_means_NGM[,vect_P_traits[k]],xlim=c(-1,2),type="n",bty="n",las=1,xaxt="n",xlab="",ylab=expression(paste("Area  (50 x ",mm^2,")")),
         main=c("SF","SB","FS","FB","BS",'BF','Size')[k])
  }

    
  #if(k>3) axis(side=1,at=c(0,1),labels=c("A6140","GA[1,2,4]50"),padj=.5)
   axis(side=1,at=c(0,1),labels=c("A6140","GA[1,2,4]50"),padj=.5)
  temp <- subset(model_line_means,population=="A6140" & env_label=='NGM')
  temp2 <- subset(model_line_means,population=="GA150" & env_label=='NGM')
  temp3 <- subset(model_line_means,population=="GA250" & env_label=='NGM')
  temp4 <- subset(model_line_means,population=="GA450" & env_label=='NGM')
  
  points(jitter(rep(posX,c(nrow(temp),nrow(temp2),nrow(temp3),nrow(temp4)))),rbind(temp,temp2,temp3,temp4)[,vect_P_traits[k]],
         pch=16,cex=.8,col=rep(c("grey",col_GA),c(nrow(temp),nrow(temp2),nrow(temp3),nrow(temp4))))
  
  
  arrows(posX , model_populations_means[model_populations_means$env_label=="NGM",paste0(vect_P_traits[k],"_lower.CL")], posX,model_populations_means[model_populations_means$env_label=="NGM",paste0(vect_P_traits[k],"_upper.CL")],code=3,angle=90,length=.1,col='black')
  points(posX, model_populations_means[model_populations_means$env_label=="NGM",paste0(vect_P_traits[k],"_emmean")],pch=21,bg="white",cex=1,col="black",type="b",lwd=2)
  
}
dev.off()
