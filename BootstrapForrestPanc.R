# This program is written to perform a bootstrap randomforrest for texture features
# to predict outcome. Version: 3.30.2021. Written by Garrett Simpson. 
# randomForest will be used to pre-select important variables using Gini Index.
# Leave-one-out analysis will be performed as follows:
#     1: select patients for inclusion
#     2: create bootstap sample of data
#     3: predict response of patient left out for prediction and record
#     4: repeat 2&3 1000 times
#     5: calculate AUC and 95% confidence intervals
#     6: repeat 1-5 for the number of patients included in study
#     7: calculate mean AUC and CI range [this AUC is called internal validation AUC]
#
# clear global environment
rm(list=ls()) 
# load relevant libraries
library('randomForest')
library('pROC')
library('psych')
library('boot')

#function to bootstrap data and test with logistic regression model
bootstrappingGLM <- function(dataIN,numBoots){
  # Function designed and written by Garrett Simpson to perform bootstrapped
  # random forest. Will perform bootstrap 'numBoots' times using 'data'.
  #Input  Data: Data column 1 should be binary response levels
  #             Data column 2 should be biniary response (0/1)
  #             column 2 and 3 should be texture features for prediction/evaluation
  #Output/return: AUC mean with 25 and 75 percentile ranges
  #             (internal validation auc)
  #             data will be split randomly into 2/3 for model training
  #             AUC/model is evaluated using all of data
  # load relevant libraries
  library('pROC')
  library('caTools')
  #testing data
  #tdata <- read.csv("C:\\Users\\gns24\\Documents\\PancE64_28pts_2.10.2021.csv",header=TRUE)
  #dataIN <-tdata[,2:5]
  
  dataIN<-data.frame(responsed=dataIN[,1],binaryresponse=dataIN[,2],feature1=dataIN[,3],feature2=dataIN[,4])
  nmbrPts<-length(dataIN[,1])
  pts<-1:nmbrPts[1]
  aucVector<-NULL
  for(bsn in 1:1000){
    sample<-NULL
    sample <- sample.split(pts, SplitRatio = 2/3)
    trainingpts<-NULL
    trainingpts<- subset(pts, sample == TRUE)
    subdata <- dataIN[trainingpts,1]
    LRModel<-NULL
    preds<-NULL
    LRModel <-glm(formula = responsed~feature1+feature2,data=dataIN[trainingpts,],family=binomial(logit))
    preds <- predict.glm(LRModel,dataIN[,3:4],type="response")
    tmpPerform<-roc(dataIN[,1],preds)
    aucVector[bsn]<-as.numeric(auc(tmpPerform))
  }
  
  return(aucVector)
  
}


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

BED40FeatureMdl <- randomForest(x=BED40,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
BED40imp <- importance(BED40FeatureMdl)
BED40orderImp <- order(BED40imp[,4],decreasing=TRUE)
BED40TopSix <- textureName[BED40orderImp[1:6]]

# Top features for each
topfeatures <- cbind(BED20TopSix,BED40TopSix)

# Testing top two feature's performance using LOO for each amount of BED

r<-c(1,2)
#Select which type of features you want to run using the leave-one-out
#reference top 2 features automatically by variable name
dataLOO20 <- BED20[c(BED20TopSix[r])]
dataLOO40 <- BED40[c(BED40TopSix[r])]
varImpPlot(BED20FeatureMdl,sort=TRUE,n.var=10,type=2)
varImpPlot(BED40FeatureMdl,sort=TRUE,n.var=10,type=2)


responses <- PancData[,3]
binaryRepsonse <- PancData[,2]

bootsamples<-1000
dataLOO1in <- cbind(responses,binaryRepsonse,dataLOO20)
dataLOO1AUC <- bootstrappingGLM(dataLOO1in,bootsamples)
FX20AUC<-mean(dataLOO1AUC)
quantiles205 <- quantile(dataLOO1AUC,probs=seq(0,1,.025), type=3,na.rm =TRUE)

dataLOO2in <- cbind(responses,binaryRepsonse,dataLOO40)
dataLOO2AUC <- bootstrappingGLM(dataLOO2in,bootsamples)
FX40AUC<-mean(dataLOO2AUC)
quantiles405 <- quantile(dataLOO2AUC,probs=seq(0,1,.025), type=3,na.rm =TRUE)

quantpositions<-c(2,40)
print(FX20AUC)
print(quantiles205[quantpositions])
print(FX40AUC)
print(quantiles405[quantpositions])



nmbrPts<-dim(PancData)
pts<-1:nmbrPts[1]
pred<-NULL
responsepred <-NULL
truerespsonse<-NULL
r<-c(1,2)
#Select which type of features you want to run using the leave-one-out
#reference top 2 features automatically by variable name
dataLOO20 <- BED20[c(BED20TopSix[r])]
dataLOO40 <- BED40[c(BED40TopSix[r])]
responses <- PancData[,3]
#testing BED features selected from liver in the pancreas data
rtemp <- c(1,23)
dataLOO40 <- BED40[rtemp]

print(cor(fx1[c(1)],fx1[c(23)],method="spearman"))
print(cor(fx1[c(1)],fx1[c(37)],method="spearman"))
print(cor(fx1[c(23)],fx1[c(38)],method="spearman"))
print(cor(fx1[c(23)],fx1[c(21)],method="spearman"))


# nmbrPts<-dim(PancData)
# pts<-1:nmbrPts[1]
# responsepred <-NULL
# respreal <-NULL
# pred<-NULL
# resptemp<-NULL
# rngLOOpt<-1:nmbrPts[1]
# for (ii in 1:1000){
#   pt_to_exclude<- sample(rngLOOpt,size = 1,replace = TRUE)
#   #print(pt_to_exclude)
#   ptstesting<- pts[pts!=pt_to_exclude]
#   #print(ptstesting)
#   ptresponses <- PancData[ptstesting,3]
#   #print(ptresponses)
#   mdlDefault<-NULL
#   mdlDefault <- randomForest(x=dataLOO[ptstesting,],y=ptresponses,ntree=500,mtry=2)
#   predf <- predict(mdlDefault,dataLOO[pt_to_exclude,],type='response')
#   #print(predf)
#   pred[ii]<- predf
#   #print(pred)
#   respreal[ii]<- PancData[pt_to_exclude,3]
#   #print(respreal)
# }
# rocobj<-roc(response=respreal,predictor=pred,ci="true")
# print(rocobj)
# plot(rocobj)
