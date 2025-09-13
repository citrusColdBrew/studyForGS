#Import data
rm(list=ls())
options(repos="http://mirrors.ustc.edu.cn/CRAN/")
setwd("/Users/sameen/workspace/statistical genomics/Lecture01_phenotype")

source("gapit_functions.txt")

myGD=read.table("mdp_numeric.txt",head=T)
#mdp_numeric.txt 是数值型基因型数据 (GD)，其中第一列是品系名称 (taxa)，其余列是SNP标记，纯合子用"0"和"2"表示，杂合子用"1"表示。head=T 表示文件的第一行是列头。
myGM=read.table("mdp_SNP_information.txt",head=T)
#mdp_SNP_information.txt 包含了基因型图谱信息 (GM)，即SNP的名称、染色体编号和物理位置信息。在GAPIT分析中，GD 和 GM 文件中的SNP必须保持一致的顺序。

#Import function
#setwd("/Users/Jiabo/Documents/China/SWUN/Lesson/animal-experiment-design/Materials")
source("G2P.R")
source("GWASbyCor.R")

#Simulation
X=myGD[,-1]#移除第一列
index1to5=myGM[,2]<6#创建一个逻辑向量 index1to5。myGM[,2] 是染色体编号列。这个向量将标记 myGM 中染色体编号小于6（即前5条染色体）的SNP为 TRUE
X1to5 = X[,index1to5]#根据 index1to5 筛选 X，只保留位于前5条染色体上的SNP数据。
set.seed(99164)#设置随机数种子。这使得随机模拟的结果是可重复的，每次运行代码都会得到相同的随机结果。
mySim=G2P(X= X1to5,h2=.75,alpha=1,NQTN=10,distribution="norm")

#GWAS
# p=GWASbyCor(X=X,y=mySim$y)

# color.vector <- rep(c("deepskyblue","orange","forestgreen","indianred3"),10)
# m=nrow(myGM)
# plot(t(-log10(p))~seq(1:m),col=color.vector[myGM[,2]])
# abline(v=mySim$QTN.position, lty = 2, lwd=2, col = "black")
myY=cbind(as.data.frame(myGD[,1]),mySim$y)


myGAPIT=GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
QTN.position=mySim$QTN.position,
PCA.total=3,
Random.model=FALSE,
file.output=T,
model=c("GLM", "MLM", "CMLM", "SUPER", "MLMM", "FarmCPU", "Blink"))


library(parallel)
library(doParallel)
source("http://zzlab.net/GAPIT/gapit_functions.txt")

# Import demo data
myGD <- read.table(file = "http://zzlab.net/GAPIT/data/mdp_numeric.txt", head = TRUE)
myGM <- read.table(file = "http://zzlab.net/GAPIT/data/mdp_SNP_information.txt", head = TRUE)

# Simulate 10 QTN on the first half chromosomes将模拟10个数量性状核苷酸 (QTN) 位点，并且这些位点将位于前一半的染色体上。QTN是影响数量性状的基因组区域。
index1to5 <- myGM[, 2] < 6
set.seed(99164)
mySim <- GAPIT.Phenotype.Simulation(GD = myGD[, c(TRUE, index1to5)], GM = myGM[index1to5, ],
                                   h2 = 0.7, NQTN = 40, effectunit = 0.95, QTNDist = "normal")

# GWAS with GAPIT using parallel processing
num_cores <- detectCores()  # Detect the number of available CPU cores
cl <- makeCluster(num_cores)  # Create a cluster with all CPU cores
models = c("GLM", "MLM", "CMLM", "SUPER", "MLMM", "FarmCPU", "Blink")
registerDoParallel(cl)  # Register the cluster for parallel processing

myGAPIT <- foreach(model = c("GLM", "MLM", "CMLM", "SUPER", "MLMM", "FarmCPU", "Blink"), .combine = c) %dopar% {
  GAPIT(Y = mySim$Y, GD = myGD, GM = myGM, PCA.total = 3,
        QTN.position = mySim$QTN.position, model = model)
}

stopCluster(cl)  # Stop the cluster


### Output multiple Manhattans with GWAS resutls, which were read from the located direction

   GMM=GAPIT.Multiple.Manhattan(model_store=models,Y.names=colnames(mySim$Y)[-1],GM=myGM,seqQTN = mySim$QTN.position,plot.type=c("w","h"))
   GAPIT.Circle.Manhattan.Plot(band=1,r=3,GMM$multip_mapP,plot.type=c("c","q"),signal.line=1,xz=GMM$xz,threshold=0.01)
   GMM=GAPIT.Multiple.Manhattan(model_store=models,Y.names=colnames(mySim$Y)[-1],GM=myGM,seqQTN = mySim$QTN.position,plot.type=c("s"))

########GS
setwd("/Users/sameen/workspace/statistical genomics/Lecture01_phenotype")

source("gapit_functions.txt")

myGD=read.table("mdp_numeric.txt",head=T)
myGM=read.table("mdp_SNP_information.txt",head=T)

#Simultate 10 QTN on the first half chromosomes
X=myGD[,-1]
index1to5=myGM[,2]<6
X1to5 = X[,index1to5]
taxa=myGD[,1]#再次提取 myGD 中的第一列，即品系名称，赋值给 taxa 变量

set.seed(99164)
GD.candidate=cbind(taxa,X1to5)#将品系名称 (taxa) 和筛选后的基因型数据 (X1to5) 横向合并，形成用于模拟的候选基因型数据 GD.candidate。

mySim<-GAPIT(h2=0.75,NQTN=20,GD= GD.candidate,GM=myGM[index1to5,],PCA.total=3)
#调用 GAPIT 主函数 进行表型模拟。
#  ◦ h2=0.75: 设置模拟的狭义遗传力为0.75。
#  ◦ NQTN=20: 设置模拟的QTN数量为20个。
#  ◦ GD= GD.candidate: 使用前面准备的候选基因型数据。
#  ◦ GM=myGM[index1to5,]: 使用筛选后的基因型图谱信息。
#  ◦ PCA.total=3: 指定在分析中包含3个主成分作为协变量。
#  ◦ 这个函数调用将模拟出一个性状，并将其结果（包括模拟的表型和QTN位置）存储在 mySim 对象中。

posi=mySim$QTN.position#从模拟结果 mySim 中提取模拟QTN的实际位置，赋值给 posi 变量。
myY=mySim$Y#从模拟结果 mySim 中提取模拟的表型数据，赋值给 myY 变量。

myGAPIT <- GAPIT(
Y=mySim$Y,
GD=myGD,
GM=myGM,
PCA.total=3,
Random.model=FALSE,
QTN.position=mySim$QTN.position,
model="GLM")
pred=merge(mySim$Y,myGAPIT$Pred,by.x="Taxa",by.y="Taxa")
u=cbind(as.data.frame(mySim$Y[,1]),mySim$u)
colnames(u)=c("Taxa","u")
pred.u=merge(u,myGAPIT$Pred,by.x="Taxa",by.y="Taxa")
ry2=cor(pred[,2],pred[,12])^2
ru2=cor(pred.u[,2],pred.u[,12])^2
par(mfrow=c(2,1), mar = c(3,4,1,1))
plot(pred[,2],pred[,12])
mtext(paste("R square=",ry2,sep=""), side = 3)
plot(pred.u[,2],pred.u[,12])
mtext(paste("R square=",ru2,sep=""), side = 3)


ntop=10
index=order(myGAPIT$GWAS$P.value)
top=index[1:ntop]
myQTN=cbind(myGAPIT$PCA[,1:4], myGD[,c(top+1)])

myGAPIT2<- GAPIT(
Y=mySim$Y,
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
ry2=cor(pred[,2],pred[,12])^2
ru2=cor(pred.u[,2],pred.u[,12])^2
par(mfrow=c(2,1), mar = c(3,4,1,1))
plot(pred[,2],pred[,12])
mtext(paste("R square=",ry2,sep=""), side = 3)
plot(pred.u[,2],pred.u[,12])
mtext(paste("R square=",ru2,sep=""), side = 3)

ntop=200
index=order(myGAPIT$GWAS$P.value)
top=index[1:ntop]
myQTN=cbind(myGAPIT$PCA[,1:4], myGD[,c(top+1)])

myGAPIT3<- GAPIT(
Y=mySim$Y,
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
ry2=cor(pred[,2],pred[,12])^2
ru2=cor(pred.u[,2],pred.u[,12])^2
par(mfrow=c(2,1), mar = c(3,4,1,1))
plot(pred[,2],pred[,12])
mtext(paste("R square=",ry2,sep=""), side = 3)
plot(pred.u[,2],pred.u[,12])
mtext(paste("R square=",ru2,sep=""), side = 3)

#####cross validation




set.seed(99163)
mysimulation<-GAPIT.Phenotype.Simulation(h2=0.75,NQTN=20,GD=myGD,GM=myGM)

posi=mysimulation$QTN.position
myY=mysimulation$Y
nfold=5

ref.acc=NULL
inf.acc=NULL
sets=sample(cut(1:nrow(myY ),nfold,labels=FALSE),nrow(myY ))

for(j in 1:nfold)
{
   training=myY[,c(1,2)]
   training[sets==j,2]=NA
   training_index=is.na(training[,2])
   testing=myY[training_index,c(1,2)]
   colnames(myY)=c("Taxa","Trait")

   myGAPIT=GAPIT(
       Y=training,
       GD=myGD,
	   GM=myGM,
   	   PCA.total=3,
	   model="gBLUP",
	   file.output=T)

   real.prediction=merge(myY,myGAPIT$Pred,by.x="Taxa",by.y="Taxa")
   inference.r=cor(real.prediction[training_index,2],real.prediction[training_index,12])
   reference.r=cor(real.prediction[!training_index,2],real.prediction[!training_index,12])

ref.acc=append(ref.acc,reference.r)
inf.acc=append(inf.acc,inference.r)
}

ref.acc

inf.acc
