rm(list = ls())
gc()
 library(MCMCglmm)

# This file export the eigendecomposition of the G matrices

load("output_files/RData/G_matrices_high_salt.RData")
load("output_files/RData/G_matrices_low_salt.RData")
vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")

### Output the matrices

ev_NaCl_A6 = NULL
ev_NaCl_GA1 = NULL
ev_NaCl_GA2 = NULL
ev_NaCl_GA4 = NULL

ev_NGM_A6 = NULL
ev_NGM_GA1 = NULL
ev_NGM_GA2 = NULL
ev_NGM_GA4 = NULL

for(i in 1:nrow(VCV_mat_NaCl[[1]]$VCV_Mat)){
  
  ev_NaCl_A6=rbind(ev_NaCl_A6,eigen(matrix(VCV_mat_NaCl[[1]]$VCV_Mat[i,1:49],7,7)/2)$values)
  ev_NaCl_GA1=rbind(ev_NaCl_GA1,eigen(matrix(VCV_mat_NaCl[[2]]$VCV_Mat[i,1:49],7,7)/2)$values)
  ev_NaCl_GA2=rbind(ev_NaCl_GA2,eigen(matrix(VCV_mat_NaCl[[3]]$VCV_Mat[i,1:49],7,7)/2)$values)
  ev_NaCl_GA4=rbind(ev_NaCl_GA4,eigen(matrix(VCV_mat_NaCl[[4]]$VCV_Mat[i,1:49],7,7)/2)$values)
  
  ev_NGM_A6=rbind(ev_NGM_A6,eigen(matrix(VCV_mat_NGM[[1]]$VCV_Mat[i,1:49],7,7)/2)$values)
  ev_NGM_GA1=rbind(ev_NGM_GA1,eigen(matrix(VCV_mat_NGM[[2]]$VCV_Mat[i,1:49],7,7)/2)$values)
  ev_NGM_GA2=rbind(ev_NGM_GA2,eigen(matrix(VCV_mat_NGM[[3]]$VCV_Mat[i,1:49],7,7)/2)$values)  
  ev_NGM_GA4=rbind(ev_NGM_GA4,eigen(matrix(VCV_mat_NGM[[4]]$VCV_Mat[i,1:49],7,7)/2)$values)    
  }

list_all_ev <- list(ev_NaCl_A6,ev_NaCl_GA1,ev_NaCl_GA2,ev_NaCl_GA4,ev_NGM_A6,ev_NGM_GA1,ev_NGM_GA2,ev_NGM_GA4)

for(i in 1:length(list_all_ev)){
  
temp=round(t(cbind(posterior.mode(as.mcmc(list_all_ev[[i]])),HPDinterval(as.mcmc(list_all_ev[[i]])))),digits=3)
temp2 = rbind(temp,round(temp[1,]/sum(temp[1,]),digits=3))

if(i%in%c(1:4)){
temp3 = as.data.frame(rbind(temp2,round(eigen(VCV_mat_NaCl[[i]]$G1_mat/2)$vectors,digits=3)))
}else{
  temp3 = as.data.frame(rbind(temp2,round(eigen(VCV_mat_NGM[[i-4]]$G1_mat/2)$vectors,digits=3)))
}
names(temp3)=paste0("g",c("max",2:7))
rownames(temp3)=c("Eigenvalues","HPD lower","HPD upper","Proportion",vect_P_traits)

file.name=paste0("output_files/G_matrix_eigendecomposition/",c("A6140","GA150","GA250",'GA450'),rep(c("_NaCl","_NGM"),each=4),'.csv')[i]
write.csv(temp3,file=file.name)

}