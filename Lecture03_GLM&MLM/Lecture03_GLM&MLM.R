rm(list=ls())
#Import data
rm(list=ls())
options(repos="http://mirrors.ustc.edu.cn/CRAN/")
source("gapit_functions.R")
#setwd("Lecture03_GLM&MLM/")

myGD=read.table("Lecture01_phenotype/mdp_numeric.txt",head=T)
myGM=read.table("Lecture01_phenotype/mdp_SNP_information.txt",head=T)

source("Lecture01_phenotype/G2P.R")
source("Lecture01_phenotype/GWASbyCor.R")

#Simulation
X=myGD[,-1]
index1to5=myGM[,2]<6
X1to5 = X[,index1to5]
set.seed(99164)
mySim=G2P(X= X1to5,h2=.75,alpha=1,NQTN=10,distribution="norm")

#GWAS
p= GWASbyCor(X=X,y=mySim$y)

color.vector <- rep(c("deepskyblue","orange","forestgreen","indianred3"),10)
m=nrow(myGM)
plot(t(-log10(p))~seq(1:m),col=color.vector[myGM[,2]])
abline(v=mySim$QTN.position, lty = 2, lwd=2, col = "black")
myY=cbind(as.data.frame(myGD[,1]),mySim$y)

p.obs=p
m2=length(p.obs)
p.uni=runif(m2,0,1)
order.obs=order(p.obs)
order.uni=order(p.uni)

plot(-log10(p.uni[order.uni]),
-log10(p.obs[order.obs]), )
abline(a = 0, b = 1, col = "red")

PCA=prcomp(X)
plot(mySim$y,PCA$x[,2])
cor(mySim$y,PCA$x[,2])


y=mySim$y
G=myGD[,-1]
n=nrow(G)
m=ncol(G)
P=matrix(NA,1,m)

for (i in 1:m){
x=G[,i]
if(max(x)==min(x)){
p=1}else{
X=cbind(1, PCA$x[,2],x)
LHS=t(X)%*%X
C=solve(LHS)
RHS=t(X)%*%y
b=C%*%RHS
yb=X%*%b
e=y-yb
n=length(y)
ve=sum(e^2)/(n-1)
vt=C*ve
t=b/sqrt(diag(vt))
p=2*(1-pt(abs(t),n-2))
} #end of testing variation
P[i]=p[length(p)]
} #end of looping for markers


p.obs=P
m2=length(p.obs)
p.uni=runif(m2,0,1)
order.obs=order(p.obs)
order.uni=order(p.uni)

plot(-log10(p.uni[order.uni]),
-log10(p.obs[order.obs]), )
abline(a = 0, b = 1, col = "red")


G=myGD[,-1]
n=nrow(G)
m=ncol(G)
P=matrix(NA,1,m)
for (i in 1:m){
x=G[,i]
if(max(x)==min(x)){
p=1}else{
X=cbind(1, PCA$x[,1:3],x)
LHS=t(X)%*%X
C=solve(LHS)
RHS=t(X)%*%y
b=C%*%RHS
yb=X%*%b
e=y-yb
n=length(y)
ve=sum(e^2)/(n-1)
vt=C*ve
t=b/sqrt(diag(vt))
p=2*(1-pt(abs(t),n-2))
} #end of testing variation
P[i]=p[length(p)]
} #end of looping for markers


p.obs=P
m2=length(p.obs)
p.uni=runif(m2,0,1)
order.obs=order(p.obs)
order.uni=order(p.uni)

plot(-log10(p.uni[order.uni]),
-log10(p.obs[order.obs]), )
abline(a = 0, b = 1, col = "red")


color.vector <- rep(c("deepskyblue","orange","forestgreen","indianred3"),10)
m=nrow(myGM)
plot(t(-log10(P))~seq(1:m),col=color.vector[myGM[,2]])
abline(v=mySim$QTN.position, lty = 2, lwd=2, col = "black")


myGAPIT=GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
QTN.position=mySim$QTN.position,
PCA.total=3,
file.output=FALSE,
model="MLM")
