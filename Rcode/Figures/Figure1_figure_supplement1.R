rm(list=ls())
gc()

coef_manova <- read.table('output_files/txt/Manova.coefficients.txt',sep="\t")
coef_manova_GA_plast = coef_manova[1,]
coef_univ_plast = -read.table(file='output_files/txt/Plasticity_contrasts.txt',sep='\t',h=TRUE)[,2]

pdf("plots/Figure1_figure_supplement1.pdf",w=5,h=5)
par(mar=c(5,5,4,2))
plot(coef_univ_plast,as.matrix(coef_manova_GA_plast),pch=16,col="black",xlab="Estimates from \nunivariate models",ylab='Estimates \nfrom MANOVA ',las=1,bty="n")
abline(a=0,b=1)
dev.off()


