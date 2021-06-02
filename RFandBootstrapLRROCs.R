rm(list=ls()) 
# load relevant libraries
library('randomForest')
library('pROC')
library('psych')
library('boot')
library('pROC')
library('caTools')

# Location of CSV with data organized as described above
PancData <- read.csv("C:\\Users\\gns24\\Documents\\PancE64_30ptsBED_4.29.2021.csv",header=TRUE)
# Number of patients
npatients <- 1:nrow(data)
# Data is already organized by BED delivered per fraction for calculation of
# Texture feature changes as a function of dose
fx1 <- PancData[,7:45]
fx2 <- PancData[,46:84]
fx3 <- PancData[,85:123]
# Feature change calculation
BED20 <- fx1-fx2
BED40 <- fx1-fx3
# Defining response variables and coding for use in randomforest algorithm
labels <- c("NR","RS")
levels <- c(0,1)
ptResponses <- factor(x=PancData[,3],levels,labels)
# Texture names
textureName<-colnames(fx1)
# Top six features for 20 and 40 bed
BED20FeatureMdl <- randomForest(x=BED20,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
BED20imp <- importance(BED20FeatureMdl)
BED20orderImp <- order(BED20imp[,4],decreasing=TRUE)
BED20TopSix <- textureName[BED20orderImp[1:6]]
#varImpPlot(BED20FeatureMdl,sort=TRUE,n.var=10,type=2)


BED40FeatureMdl <- randomForest(x=BED40,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
BED40imp <- importance(BED40FeatureMdl)
BED40orderImp <- order(BED40imp[,4],decreasing=TRUE)
BED40TopSix <- textureName[BED40orderImp[1:6]]
#varImpPlot(BED40FeatureMdl,sort=TRUE,n.var=10,type=2)
# Top features for each
topfeatures <- cbind(BED20TopSix,BED40TopSix)

# Testing top two feature's performance using LOO for each amount of BED

r<-c(1,2)
#Select which type of features you want to run using the leave-one-out
#reference top 2 features automatically by variable name
dataLOO20 <- BED20[c(BED20TopSix[r])]
dataLOO40 <- BED40[c(BED40TopSix[r])]
responses <- PancData[,3]
binaryRepsonse <- PancData[,2]
bootsamples<-1000
nmbrPts<-dim(PancData)
pts<-1:nmbrPts[1]
pred<-NULL
responsepred <-NULL
truerespsonse<-NULL
AUC20BED<-NULL
AUC40BED<-NULL

pdf("BED20.pdf")
dataLOO20 <- BED20[c(BED20TopSix[r])]
dataLOO40 <- BED40[c(BED40TopSix[r])]
bed20mdlroc <- randomForest(x=dataLOO20,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=2,header=TRUE)
bed40mdlroc<- randomForest(x=dataLOO40,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=2,header=TRUE)
preds20 <-as.numeric(predict(bed20mdlroc,dataLOO20, type = "response"))
preds40 <-as.numeric(predict(bed40mdlroc,dataLOO40, type = "response"))
aucVector20<-NULL
aucVector40<-NULL
rocobj<-roc(response=PancData[,2],predictor=bed20mdlroc$votes[,2])

plot.roc(rocobj)
dataIN<-NULL
dataIN<-data.frame(responsed=PancData[,2],binaryresponse=PancData[,3],feature1=dataLOO20[,1],feature2=dataLOO20[,2])
check<-NULL
count<-0
for (modl20 in 1:bootsamples){
  sample<-NULL
  sample <- sample.split(pts, SplitRatio = 2/3)
  trainingpts<-NULL
  trainingpts<- subset(pts, sample == TRUE)
  LRModel<-NULL
  preds<-NULL
  tmpPerform<-NULL
  LRModel <-glm(formula = responsed~feature1+feature2,data=dataIN[trainingpts,],family=binomial(logit), maxit = 100)
  preds <- predict.glm(LRModel,dataIN[,3:4],type="response")
  tmpPerform<-roc(dataIN[,1],preds)
  plot.roc(tmpPerform,add=TRUE,col = '#D3D3D3')
}
plot.roc(rocobj,add=TRUE)
dev.off()

pdf("BED40.pdf")
rocobj40<-roc(response=PancData[,2],predictor=bed40mdlroc$votes[,2])
plot(rocobj40)
dataIN<-NULL
dataIN<-data.frame(responsed=PancData[,2],binaryresponse=PancData[,3],feature1=dataLOO40[,1],feature2=dataLOO40[,2])

for (modl40 in 1:bootsamples){
  sample<-NULL
  sample <- sample.split(pts, SplitRatio = 2/3)
  trainingpts<-NULL
  trainingpts<- subset(pts, sample == TRUE)
  subdata <- dataIN[trainingpts,1]
  LRModel<-NULL
  preds<-NULL
  tmpPerform<-NULL
  LRModel <-glm(formula = responsed~feature1+feature2,data=dataIN[trainingpts,],family=binomial(logit))
  preds <- predict.glm(LRModel,dataIN[,3:4],type="response")
  tmpPerform<-roc(dataIN[,1],preds)
  aucVector40[modl40]<-as.numeric(auc(tmpPerform))
  plot.roc(tmpPerform,add=TRUE,col = '#D3D3D3')
}
plot(rocobj40,add=TRUE)
dev.off()

