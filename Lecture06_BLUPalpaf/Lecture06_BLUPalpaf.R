
setwd("/Users/sameen/workspace/statistical genomics/code_day2")
# reading data and source functions
source("gapit_functions.txt")
myGD=read.table(file="mdp_numeric.txt",head=T) 
myGM=read.table(file="mdp_SNP_information.txt",head=T)
# divide population for ref and inf
set.seed(99164)
n=nrow(myGD)
testing=sample(n,round(n/5),replace=F)
training=-testing

#Simultate 20 QTN on the first half chromosomes
X=myGD[,-1]
taxa=myGD[,1]

set.seed(99164)

mySim<-GAPIT(h2=0.75,NQTN=20,GD= myGD,GM=myGM)
posi=mySim$QTN.position
myY=mySim$Y


myGAPIT5 <- GAPIT(
Y=mySim$Y[training,],
GD=myGD, 
GM=myGM,
PCA.total=3,
model="gBLUP", 
SNP.test=FALSE)


pred=merge(mySim$Y,myGAPIT5$Pred,by.x="Taxa",by.y="Taxa")
u=cbind(as.data.frame(mySim$Y[,1]),mySim$u)
colnames(u)=c("Taxa","u")
pred.u=merge(u,myGAPIT5$Pred,by.x="Taxa",by.y="Taxa")
ry2=cor(pred[testing,2],pred[testing,12])^2
ru2=cor(pred.u[testing,2],pred.u[testing,12])^2
par(mfrow=c(2,1), mar = c(3,4,1,1))
plot(pred[testing,2],pred[testing,12])
mtext(paste("R square=",ry2,sep=""), side = 3)
plot(pred.u[testing,2],pred.u[testing,12])
mtext(paste("R square=",ru2,sep=""), side = 3)


myGAPIT6 <- GAPIT(
Y=mySim$Y[training,],
GD=myGD, 
GM=myGM,
PCA.total=3,
model="cBLUP", 
SNP.test=FALSE)


pred=merge(mySim$Y,myGAPIT6$Pred,by.x="Taxa",by.y="Taxa")
u=cbind(as.data.frame(mySim$Y[,1]),mySim$u)
colnames(u)=c("Taxa","u")
pred.u=merge(u,myGAPIT6$Pred,by.x="Taxa",by.y="Taxa")
ry2=cor(pred[testing,2],pred[testing,12])^2
ru2=cor(pred.u[testing,2],pred.u[testing,12])^2
par(mfrow=c(2,1), mar = c(3,4,1,1))
plot(pred[testing,2],pred[testing,12])
mtext(paste("R square=",ry2,sep=""), side = 3)
plot(pred.u[testing,2],pred.u[testing,12])
mtext(paste("R square=",ru2,sep=""), side = 3)


myGAPIT7 <- GAPIT(
Y=mySim$Y[training,],
GD=myGD, 
GM=myGM,
PCA.total=3,
model="sBLUP", 
SNP.test=FALSE)


pred=merge(mySim$Y,myGAPIT7$Pred,by.x="Taxa",by.y="Taxa")
u=cbind(as.data.frame(mySim$Y[,1]),mySim$u)
colnames(u)=c("Taxa","u")
pred.u=merge(u,myGAPIT7$Pred,by.x="Taxa",by.y="Taxa")
ry2=cor(pred[testing,2],pred[testing,12])^2
ru2=cor(pred.u[testing,2],pred.u[testing,12])^2
par(mfrow=c(2,1), mar = c(3,4,1,1))
plot(pred[testing,2],pred[testing,12])
mtext(paste("R square=",ry2,sep=""), side = 3)
plot(pred.u[testing,2],pred.u[testing,12])
mtext(paste("R square=",ru2,sep=""), side = 3)



myGAPIT8 <- GAPIT(
Y=mySim$Y[training,],
GD=myGD, 
GM=myGM,
PCA.total=3,
model="BLINK", 
buspred=TRUE,
lmpred=c(FALSE)
)


pred=merge(mySim$Y,myGAPIT8$Pred,by.x="Taxa",by.y="Taxa")
u=cbind(as.data.frame(mySim$Y[,1]),mySim$u)
colnames(u)=c("Taxa","u")
pred.u=merge(u,myGAPIT8$Pred,by.x="Taxa",by.y="Taxa")
ry2=cor(pred[testing,2],pred[testing,12])^2
ru2=cor(pred.u[testing,2],pred.u[testing,12])^2
par(mfrow=c(2,1), mar = c(3,4,1,1))
plot(pred[testing,2],pred[testing,12])
mtext(paste("R square=",ry2,sep=""), side = 3)
plot(pred.u[testing,2],pred.u[testing,12])
mtext(paste("R square=",ru2,sep=""), side = 3)
