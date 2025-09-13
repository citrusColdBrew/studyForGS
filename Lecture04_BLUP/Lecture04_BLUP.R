rm(list=ls())
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
# GD.candidate=cbind(taxa,X1to5)

mySim<-GAPIT(h2=0.75,NQTN=20,GD= myGD,GM=myGM)
posi=mySim$QTN.position
myY=mySim$Y

# GWAS
myGAPIT <- GAPIT(
Y=mySim$Y[training,],
GD=myGD,
GM=myGM,
PCA.total=3,
Random.model=FALSE,
QTN.position=mySim$QTN.position,
model="GLM")
dim(myGAPIT$Pred)

pred=merge(mySim$Y,myGAPIT$Pred,by.x="Taxa",by.y="Taxa")
u=cbind(as.data.frame(mySim$Y[,1]),mySim$u)
colnames(u)=c("Taxa","u")
pred.u=merge(u,myGAPIT$Pred,by.x="Taxa",by.y="Taxa")
ry2=cor(pred[testing,2],pred[testing,12])^2
ru2=cor(pred.u[testing,2],pred.u[testing,12])^2
par(mfrow=c(2,1), mar = c(3,4,1,1))
plot(pred[testing,2],pred[testing,12])
mtext(paste("R square=",ry2,sep=""), side = 3)
plot(pred.u[testing,2],pred.u[testing,12])
mtext(paste("R square=",ru2,sep=""), side = 3)

# use 10 top markers by GWAS

ntop=10
index=order(myGAPIT$GWAS$P.value)
top=index[1:ntop]
myQTN=cbind(myGAPIT$PCA[,1:4], myGD[,c(top+1)])

myGAPIT2<- GAPIT(
Y=mySim$Y[training,],
GD=myGD,
GM=myGM,
Random.model=FALSE,
CV=myQTN,
QTN.position=mySim$QTN.position,
model="GLM"
)

pred=merge(mySim$Y,myGAPIT2$Pred,by.x="Taxa",by.y="Taxa")
u=cbind(as.data.frame(mySim$Y[,1]),mySim$u)
colnames(u)=c("Taxa","u")
pred.u=merge(u,myGAPIT2$Pred,by.x="Taxa",by.y="Taxa")
ry2=cor(pred[testing,2],pred[testing,12])^2
ru2=cor(pred.u[testing,2],pred.u[testing,12])^2
par(mfrow=c(2,1), mar = c(3,4,1,1))
plot(pred[testing,2],pred[testing,12])
mtext(paste("R square=",ry2,sep=""), side = 3)
plot(pred.u[testing,2],pred.u[testing,12])
mtext(paste("R square=",ru2,sep=""), side = 3)

# use 200 top markers by GWAS
ntop=200
index=order(myGAPIT$GWAS$P.value)
top=index[1:ntop]
myQTN=cbind(myGAPIT$PCA[,1:4], myGD[,c(top+1)])

myGAPIT3<- GAPIT(
Y=mySim$Y[training,],
GD=myGD,
GM=myGM,
Random.model=FALSE,
CV=myQTN,
QTN.position=mySim$QTN.position,
model="GLM"
)

pred=merge(mySim$Y,myGAPIT3$Pred,by.x="Taxa",by.y="Taxa")
u=cbind(as.data.frame(mySim$Y[,1]),mySim$u)
colnames(u)=c("Taxa","u")
pred.u=merge(u,myGAPIT3$Pred,by.x="Taxa",by.y="Taxa")
ry2=cor(pred[testing,2],pred[testing,12])^2
ru2=cor(pred.u[testing,2],pred.u[testing,12])^2
par(mfrow=c(2,1), mar = c(3,4,1,1))
plot(pred[testing,2],pred[testing,12])
mtext(paste("R square=",ry2,sep=""), side = 3)
plot(pred.u[testing,2],pred.u[testing,12])
mtext(paste("R square=",ru2,sep=""), side = 3)

#GBLUP

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

library(rrBLUP)

#prepare data
y <- mySim$Y[,2]
M=as.matrix(X[,])

#Ridge Regression
ans1 <- mixed.solve(y=y,Z=M)

#gBLUP
K <- tcrossprod(M) #K = MM'
ans2 <- mixed.solve(y=y,K=K)

#Compare GEBV
plot(M%*%ans1$u, ans2$u) 



myGAPIT5 <- GAPIT(
Y=mySim$Y,
GD=myGD, 
GM=myGM,
# PCA.total=3,
model="gBLUP", 
SNP.test=FALSE)


pred=merge(mySim$Y,myGAPIT5$Pred,by.x="Taxa",by.y="Taxa")

order.raw=match(taxa,pred[,1])
pred2=pred[order.raw,]
plot(ans2$u, pred2[,9])


