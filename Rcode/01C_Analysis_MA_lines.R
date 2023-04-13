rm(list=ls());gc()
library(MCMCglmm)
library(matrixStats)
final_merged=read.table('data/MA_lines_phenotypes.txt',sep='\t',header=TRUE)
dim(final_merged) # 226 lines

vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")

### Finally we scale the body area to match the sd of the transition rates
### currently in mm^2

mean(colSds(as.matrix(final_merged[,vect_P_traits[1:6]]))) # 0.59
tapply(final_merged$area.F,final_merged$env_label,sd) # 0.01 in A6140 vs 0.0068 in MA lines

# We apply the same factor as for the inbred lines from experimental evolution
final_merged$area.F = final_merged$area.F * 50
final_merged$anc = substring(final_merged$pop_label,1,2)
table(final_merged$anc)
  
plot.acfs <- function(x) {
  n <- dim(x)[2]
  par(mfrow = c(ceiling(n/2), 2), mar = c(3, 2, 3, 0.5))
  for (i in 1:n) {
    acf(x[, i], lag.max = 100, main = colnames(x)[i], ci.type = "ma", ylim = c(-0.15, 0.15))
    grid()
  }
}

VCV_mat_MAlines = NULL
nb_trait = length(vect_P_traits)

### Normalize the env. factors
  for(j in c('temperature',"rel_humidity","logD")){
    final_merged[,j] <- (final_merged[,j]-mean(final_merged[,j]))/sd(final_merged[,j])
  }
  
  
  phen.var = diag(nb_trait) * diag(var(subset(final_merged, select = vect_P_traits)))
  prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
                    R = list(V = phen.var/3, n = nb_trait))
  
  model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32, area.F)) ~  (temperature+rel_humidity+logD)^3 + anc + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                         rcov = ~us(trait):units, 
                         family = rep("gaussian", nb_trait), data = final_merged, prior = prior_mod, verbose = FALSE,nitt=1500000, burnin=500000,thin=100)
  
  pdf(file = paste0("output_files/auto_corr_plots_MCMCglmm_MAlines/Model_pdf_MCMC.pdf"), height = 20)
  par(mfrow = c(10, 2), mar = c(2, 2, 1, 0))
  plot(model_MCMC$Sol, auto.layout = F)
  dev.off()
  
  pdf(file = paste0("output_files/auto_corr_plots_MCMCglmm_MAlines/Model_pdf_MCMC_autocorr_.pdf"), height = 10)
  plot.acfs(model_MCMC$Sol)
  dev.off()
  
  post_dist = posterior.mode(model_MCMC$VCV)
  VCV_mat_MAlines=list(Population = "MA lines", N_measurement = nrow(final_merged), G1_mat = matrix(post_dist[1:nb_trait^2], 
                                                                                      nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
                    R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)

###Create tables
  write.table(matrix(paste0(round(1000*VCV_mat_MAlines$G1_mat)/1000
                            ,"  [",round(1000*matrix(HPDinterval(VCV_mat_MAlines$VCV_Mat)[1:49,1],ncol=7))/1000,";",
                            round(1000*matrix(HPDinterval(VCV_mat_MAlines$VCV_Mat)[1:49,2],ncol=7))/1000,"]"),ncol=7),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,file=paste0("output_files/G_mat_tables/MAlines.txt"))
  
  write.table(VCV_mat_MAlines$G1_mat, quote=FALSE,sep="\t", row.names=FALSE,col.names=FALSE, file=paste0("output_files/txt/", "MAlines",".txt"))
  
## Eigendecomposition

  # Proportion of variance explained by each axis
  eigen(VCV_mat_MAlines$G1_mat)$values/sum(eigen(VCV_mat_MAlines$G1_mat)$values)
  write.table(eigen(VCV_mat_MAlines$G1_mat)$vectors,"output_files/txt/MAlines_eigenvectors.txt",row.names=FALSE,quote=FALSE,sep='\t')
  
  
  write.table(round(eigen(VCV_mat_MAlines$G1_mat)$vectors,digits=2),"tempMA.txt",row.names=FALSE,quote=FALSE,sep='\t')