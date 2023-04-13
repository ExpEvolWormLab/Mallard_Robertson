rm(list=ls());gc()

load('output_files/RData/Gqw_High_Salt.RData')
load('output_files/RData/G_matrices_high_salt.RData')

nb_trait=8

pdf(file='plots_elife/Figure4_supplement_figure1.pdf',h=9,w=5.5)

par(mar=c(5,7,4,2))

vect_6t_v2 <- c(c(2:7,11:15,20:23,29:31,38,39,47),c(57:63),c(9*(1:8)-8))
vectX_6t_v2 <- c(VCV_with_w_NaCl$G1_mat/2)[vect_6t_v2]

vProb <- .95

vect_y_v2 =c(c(11:31),c(1:7),c(35:41),8)
vect_col_v2=(c(rep("lightblue",21),rep("darkgreen",7),rep("lightblue",7),"darkgreen"))

plot(vect_6t_v2,vect_y_v2,ylim=c(44,10),yaxt="n",bty="n",xlim=c(-.4,.7),xlab="Genetic (co-)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)

lines(c(0,0),c(31.5,10.5))
lines(c(0,0),c(34.5,41.5),col="black")

axis(side=1,pos=44)

axis(side=2,at=c(11:31,35:41),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF","SF*Area",
                                           "SB*FS","SB*FB","SB*BS","SB*BF","SB*Area",
                                           "FS*FB","FS*BS","FS*BF","FS*Area",
                                           "FB*BS","FB*BF","FB*Area","BS*BF","BS*Area","BF*Area",
                                           "SF","SB","FS","FB","BS","BF","Area"),las=1)

temp_95 <- HPDinterval(VCV_with_w_NaCl$VCV_Mat[,1:nb_trait^2]/2,prob=.95)
temp_80 <- HPDinterval(VCV_with_w_NaCl$VCV_Mat[,1:nb_trait^2]/2,prob=.8)

arrows(temp_95[vect_6t_v2,1],vect_y_v2,temp_95[vect_6t_v2,2],vect_y_v2,code=3,length=.05,angle=90,lwd=2)
arrows(temp_80[vect_6t_v2,1],vect_y_v2,temp_80[vect_6t_v2,2],vect_y_v2,code=3,length=0,angle=90,lwd=3,col="grey")
points(vectX_6t_v2,vect_y_v2,pch=21,bg="black")

### Add the A6140 G matrix

vect_Var <- c(2:7,10:14,18:21,26:28,34,35,42,1,9,17,25,33,41,49)
vProb <- .95
vectX_NaCl <- c(VCV_mat_NaCl[[1]]$G1_mat/2)[vect_Var]
true_G_NaCl=VCV_mat_NaCl[[1]]$G1_mat

temp_95 <- HPDinterval(VCV_mat_NaCl[[1]]$VCV_Mat[,1:49]/2,prob=.95)
arrows(temp_95[vect_Var,1],c(11:31,35:41)+.3,temp_95[vect_Var,2],c(11:31,35:41)+.3,code=3,length=.03,angle=90,lwd=2)

temp_80 <- HPDinterval(VCV_mat_NaCl[[1]]$VCV_Mat[,1:49]/2,prob=.83)
arrows(temp_80[vect_Var,1],c(11:31,35:41)+.3, temp_80[vect_Var,2],c(11:31,35:41)+.3,code=3,length=0,angle=90,lwd=3,col="grey")
points(vectX_NaCl,c(11:31,35:41)+.3,pch=21,bg="white")

legend(.3,30,c(expression(G[qw]),"G"),lwd=1,pch=21,pt.bg=c("black","white"))

dev.off()
