p=seq(0, 1, .01)
n=4
k=n

pyp=dbinom(k,n,p)

theMax=pyp==max(pyp)
pMax=p[theMax]
plot(p,pyp,type="b",main=paste("Data=", pMax,sep=""))


ps=p*10 - 5 # Minus 5 is middle, 1 for the far left and 9 for far right
pd=dnorm(ps)

theMax=pd==max(pd)
pMax=p[theMax]
plot(p,pd,type="b",main=paste("Prior=", pMax,sep=""))


ppy=pd*pyp

theMax=ppy==max(ppy)
pMax=p[theMax]
plot(p,ppy,type="b",main=paste("Optimum=", pMax,sep=""))


####gibbs

gibbs=function (n=10000, r=.75, sd=1)
{
mat=matrix(ncol = 2, nrow = n)
x=-50
y=50
mat[1, ]=c(x, y)
for (i in 2:n) {
x=rnorm(n=1, mean=r*y, sd=1)
y=rnorm(n=1, mean=r*x, sd=1)
mat[i, ]=c(x, y)
}
mat
}

n= 10000
bvn<-gibbs(n,.75, sd=1)
cor(bvn)

batch=5000
ndisp=1000
xlim=c(min(bvn[,1]),max(bvn[,1]))
ylim=c(min(bvn[,2]),max(bvn[,2]))
for(i in 1:n){
  if(i==1)plot(bvn[i,1],bvn[i,2],xlim=xlim,ylim=ylim,pch=20,col="red")
  if(i<ndisp&i>1)points(bvn[i,1],bvn[i,2],pch=20)
  if(i>ndisp)points(bvn[i,1],bvn[i,2],col=floor(i/batch)+1)
  
  if(i<ndisp)Sys.sleep(1/i) 
  if(i==ndisp)Sys.sleep(2) 
  if(floor(i/batch)*batch==i) Sys.sleep(1) 
}

mySample=sample(100:10000,9900)
bvn2=bvn[mySample,]

par(mfrow=c(3,1), mar = c(3,4,1,1))
plot(bvn2[1:1000,1],type="b")
lines(bvn2[1:1000,2],col="red")

plot(bvn2[100:1000,1],type="b")
lines(bvn2[100:1000,2],col="red")

plot(bvn2[801:1000,1],type="b")
lines(bvn2[801:1000,2],col="red")


n= 1000000
bvn<-gibbs(n,.75, sd=1)

block=10
s=n/block
r=matrix(0,s,1)
for (i in 1:s) {
fregment=seq( ((i-1)*block+1),i*block)
x= bvn[fregment,]
r[i,]=cor(x[,1],x[,2])
}
par(mfrow=c(1,2), mar = c(3,4,1,1))
plot(r)
hist(r)


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
index1to5=myGM[,2]<6
X1to5 = X[,index1to5]
taxa=myGD[,1]

set.seed(99164)
GD.candidate=cbind(taxa,X1to5)

mySim<-GAPIT(h2=0.75,NQTN=20,GD= GD.candidate,GM=myGM[index1to5,],PCA.total=3)
posi=mySim$QTN.position
myY=mySim$Y


#Gene mapping (GWAS)
myGAPIT3 <- GAPIT(
  Y=mySim$Y[training,],
  GD=myGD,GM=myGM, PCA.total=3,
  model="BLINK", 
  QTN.position=mySim$QTN.position,
  memo="GWAS")

myCV=myGAPIT3$PCA



index=myGAPIT3$GWAS[,4]<0.05/length(myGAPIT3$GWAS[,4]) 
QTN.order=order(abs(mySim$effect),decreasing =T)
index= mySim$QTN.position[QTN.order]+1
myQTN=cbind(myGAPIT3$PCA,myGD[,c(FALSE,index[1:15])])

#MAS
myGAPIT4 <- GAPIT(
  Y=mySim$Y[training,],
  CV=myQTN, 
  model="GLM", 
  SNP.test=FALSE,
  memo="MAS")

order=match(mySim$Y[,1],myGAPIT4$Pred[,1])
myPred=myGAPIT4$Pred[order,]
ru2=cor(myPred[testing,10],mySim$u[testing])^2
plot(myPred[testing,10],mySim$u[testing])
mtext(paste("R square=",ru2,sep=""), side = 3)

#贝叶斯回归-A
library(BLR) 
nIter=2000              #### number of iteration
burnIn=1500             #### burnin a part of iteration
set.seed(99164)
myBLR =BLR(y=as.matrix(mySim$Y[training,2]),
           	XF=myCV[training,-1],
 	XR=as.matrix(myGD[training,-1]), 
	nIter=nIter,
           	burnIn=burnIn)

pred.inf=as.matrix(myGD[testing,-1])%*%myBLR$bR
ru2 <- cor(mySim$u[testing],pred.inf)^2
plot(mySim$u[testing],pred.inf)
mtext(paste("R square=",ru2,sep=""), side = 3)




library(BLR) 
nIter=2000              #### number of iteration
burnIn=1500             #### burnin a part of iteration
set.seed(99164)
myBLR =BLR(y=as.matrix(mySim$Y[training,2]),
           	XF=myCV[training,-1],
 	XL=as.matrix(myGD[training,-1]), 
	nIter=nIter,
           	burnIn=burnIn)

pred.inf=as.matrix(myGD[testing,-1])%*%myBLR$bL
ru2 <- cor(mySim$u[testing],pred.inf)^2
plot(mySim$u[testing],pred.inf)
mtext(paste("R square=",ru2,sep=""), side = 3)


#install.packages("BGLR") 
library(BGLR) 
nIter=2000              #### number of iteration
burnIn=1500             #### burnin a part of iteration
set.seed(99164)
myBGLR =BGLR(y=as.matrix(mySim$Y[training,2]),
             ETA=list(list(X=myCV[training,-1],model='FIXED', saveEffects=TRUE),
                      list(X=as.matrix(myGD[training,-1]),model='BRR', saveEffects=TRUE)
                      #list(X=as.matrix(myGD[training,-1]),model='BL', saveEffects=TRUE)
                      #list(X=as.matrix(myGD[training,-1]),model='BayesA', saveEffects=TRUE)
                      #list(X=as.matrix(myGD[training,-1]),model='BayesB', saveEffects=TRUE)
             ),
             nIter=nIter,
             burnIn=burnIn)

pred.inf=as.matrix(myGD[testing,-1])%*%myBGLR$ETA[[2]]$b
ru2 <- cor(mySim$u[testing],pred.inf)^2
plot(mySim$u[testing],pred.inf)
mtext(paste("R square=",ru2,sep=""), side = 3)


