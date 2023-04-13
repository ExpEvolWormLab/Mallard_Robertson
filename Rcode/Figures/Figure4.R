rm(list=ls());gc()

load('output_files/RData/Gqw_Low_Salt.RData')
load('output_files/RData/Gqw_High_Salt.RData')

vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32","area.F")


pdf("plots/Figure4.pdf",w=6)

par(mar=c(5,7,4,2))

vect_low <- c(57:63)
vect_high <- c(57:64)

vectX_high <- c(VCV_with_w_NaCl$G1_mat/2)[vect_high]
vectX_low  <- c(VCV_with_w_NGM$G1_mat/2)[vect_low]

vectX = c(vectX_high,vectX_low)

vProb <- .95
vectY =c(1:8-.2,1:7+.2)

plot(vectX,vectY,ylim=c(8,0),yaxt="n",bty="n",xlim=c(-.15,.3),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)
lines(c(0,0),c(.5,7.5),col="black")
axis(side=1,pos=8.5)
axis(side=2,at=c(1:8),labels=c("SF*w","SB*w","FS*w","FB*w","BS*w","BF*w","Size*w","w"),las=1)

# High Salt
temp_95_high <- HPDinterval(as.mcmc(VCV_with_w_NaCl$VCV_Mat)/2,prob=.95)
arrows(temp_95_high[vect_high,1],vectY[1:8],temp_95_high[vect_high,2],vectY[1:8],code=3,length=.05,angle=90)
temp_80_high <- HPDinterval(as.mcmc(VCV_with_w_NaCl$VCV_Mat)/2,prob=.83)
arrows(temp_80_high[vect_high,1],vectY[1:8],temp_80_high[vect_high,2],vectY[1:8],code=3,length=0,angle=90,lwd=2,col="gray")


# Low Salt
temp_95_low <- HPDinterval(as.mcmc(VCV_with_w_NGM$VCV_Mat)/2,prob=.95)
arrows(temp_95_low[vect_low,1],vectY[9:15],temp_95_low[vect_low,2],vectY[9:15],code=3,length=.05,angle=90)
temp_80_low <- HPDinterval(as.mcmc(VCV_with_w_NGM$VCV_Mat)/2,prob=.83)
arrows(temp_80_low[vect_low,1],vectY[9:15],temp_80_low[vect_low,2],vectY[9:15],code=3,length=0,angle=90,lwd=2,col="firebrick3")

points(vectX,vectY,bg="black",pch=16)
legend(.1,5,c("Low Salt","High Salt"),lwd=2,col=c("firebrick3","grey"))
dev.off()

## Export the data in a table :

output_df = data.frame(cbind(vectX,rbind(temp_80_high[vect_high,],temp_80_low[vect_low,]),rbind(temp_95_high[vect_high,],temp_95_low[vect_low,])))
rownames(output_df)=NULL
output_df$trait=c(c("SF*w","SB*w","FS*w","FB*w","BS*w","BF*w","Size*w","w"),c("SF*w","SB*w","FS*w","FB*w","BS*w","BF*w","Size*w"))
output_df$environment=rep(c("High Salt","Low Salt"),c(8,7))
output_df=output_df[,c(6:7,1:5)]
names(output_df)[3:7]=c("Post.mode","lower_83","upper_83","lower_95","upper_95")
output_df[,3:7]=round(output_df[,3:7],digits=3)
write.table(output_df,file=("output_files/txt/Selection_differentials.txt"),sep='\t',quote=FALSE,row.names=FALSE)
