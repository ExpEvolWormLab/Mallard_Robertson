rm(list=ls());gc()

coef_manova <- read.table('output_files/txt/Manova.coefficients.txt',sep="\t")

coef_manova_GA_evolution <- coef_manova[2:4,]
coef_univ = read.table(file='output_files/txt/Divergence_contrasts_High_Salt.txt',sep='\t',h=TRUE)
coef_univ_GA = - rbind(coef_univ[c(1,4,7,10,13,16,19),2] , coef_univ[c(1,4,7,10,13,16,19)+1,2] , coef_univ[c(1,4,7,10,13,16,19)+2,2])

pdf("plots/Figure6_figure_supplement1.pdf")
plot(coef_univ_GA,as.matrix(coef_manova_GA_evolution),pch=16,col="black",xlab="Estimates from \nunivariate models",ylab='Estimates \nfrom MANOVA ',las=1,bty="n")
abline(a=0,b=1)
dev.off()





