rm(list=ls())
gc()
library(lme4)
library(boot)
library(emmeans)
library(multcomp)
data = read.table("data/competition.txt",sep='\t',h=TRUE)#[,c(2:6)]

# GFP and non GFP have been counted in liquid and sometimes in two drops
# We add the numbers here
data$nonGFP[!is.na(data[,4])] = data$nonGFP[!is.na(data[,4])] + data[!is.na(data[,4]),4]
data$GFP[!is.na(data[,5])] = data$GFP[!is.na(data[,5])] + data[!is.na(data[,5]),5]

test_mod <- glm(cbind(data$GFP,data$nonGFP) ~ population ,family='quasibinomial',data=data); summary(test_mod)
test_mod_nopop <- glm(cbind(data$GFP,data$nonGFP) ~ 1 ,family='quasibinomial',data=data); summary(test_mod_nopop)

anova(test_mod,test_mod_nopop,test="F")
data$population=as.factor(data$population)

test_mod <- glm(cbind(data$GFP,data$nonGFP) ~ population ,family='quasibinomial',data=data); summary(test_mod)
em_competition<- emmeans(test_mod,~population)
pairs(em_competition,type='response')
summary(glht(test_mod,mcp(population="Tukey")))
# Both packages give similar results

write.table(capture.output(summary(glht(test_mod,mcp(population="Tukey")))),file='output_files/txt/Competition_contrasts.txt',row.names=FALSE,quote=FALSE,sep='\t')
col_GA <- c("cadetblue1", "cornflowerblue", "slateblue2")
vcol=c("black", col_GA)
for(i in 1:4) vcol[i]=rgb(t(col2rgb(vcol[i])/255),alpha=.8)

test_mod_noI <- glm(cbind(data$GFP,data$nonGFP) ~ population-1 ,family='quasibinomial',data=data); summary(test_mod)
data_freq <- data$GFP/(data$GFP+data$nonGFP)

pdf(file='plots_elife/Figure5_noasterisks.pdf',h=4,w=5)
par(mfrow=c(1,1),mar=c(5,7,4,0))
plot((1-data_freq)~jitter(as.numeric(as.factor(data$population))),pch=16,xlim=c(0,5),ylim=c(.8,1),col = vcol[as.numeric(as.factor(data$population))],ylab="Relative fitness to GFP tester\n(wild-type frequency)",xlab="Population",bty="n",
     xaxt="n")
axis(side=1,at=1:4,c('A6140',"GA150","GA250","GA450"))
arrows(1:4-.3,(1-inv.logit(confint(test_mod_noI))[,1]),1:4-.3,(1-inv.logit(confint(test_mod_noI))[,2]),code=3,length=.05,col="gray",angle=90)
points(1:4-.3,(1-inv.logit(coef(test_mod_noI))),bg=c("black", col_GA),col="gray",pch=21)
dev.off()

