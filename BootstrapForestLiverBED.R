# This program is written to perform a bootstrap randomforrest for texture features
# to predict local control. Version: 4.06.2021. Written by Garrett Simpson. 
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
#Getting the data for analysis
Liverdata <- read.csv("C:\\Users\\gns24\\Documents\\LiverTA_64Equal_BED.csv",na.strings = c("NA"," ",""),header=TRUE)
npatients <- 1:nrow(Liverdata)
#Bootstrapping function for logistic regression to evaluate top texture features
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


# organizing data for manipulation
# fx1data is the raw value of texture features from the csv file
fx1 <- Liverdata[,5:43]
fx2 <- Liverdata[,44:82]
fx3 <- Liverdata[,83:121]
fx4 <- Liverdata[,122:160]
fx5 <- Liverdata[,161:199]
# defining response variables and coding for use in randomforest algorithm
labels <- c("Yes","No")
levels <- c(1,2)
ptResponsesliver <- factor(x=Liverdata[,4],levels,labels)
#ptResponsesliver <- Liverdata[,3]
#Texture names
temp<-colnames(fx1)
#BED/fraction
BED<- Liverdata[,2]

#change as function of delivered dose
fx20BED <- fx1 -fx2 
fx40BED <- fx1 -fx3



BED20 <- randomForest(x=fx20BED,y=ptResponsesliver,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
BED20Imp <- importance(BED20)
BED20Ordered <- order(BED20Imp[,4],decreasing=TRUE)
BED20Topsix <- temp[BED20Ordered[1:6]]

BED40 <- randomForest(x=fx40BED,y=ptResponsesliver,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
BED40Imp <- importance(BED40)
BED40Ordered <- order(BED40Imp[,4],decreasing=TRUE)
BED40Topsix <- temp[BED40Ordered[1:6]]


# Testing top two feature's performance using LOO for each amount of BED
r<-c(1,2)
#Testing top BED features from Panc
#r<-c(6,37)
#Select which type of features you want to run using the leave-one-out
#reference top 2 features automatically by variable name
dataLOO20 <- fx20BED[c(BED20Topsix[r])]
dataLOO40 <- fx40BED[c(BED40Topsix[r])]
#Testing top BED features from panc
#dataLOO20 <- fx20BED[c(temp[r])]
#dataLOO40 <- fx40BED[c(temp[r])]

responses <- Liverdata[,3]
binaryRepsonse <- Liverdata[,4]-1

bootsamples<-1000
quantilspostitions<-c(2,40)
dataLOO1in <- cbind(responses,binaryRepsonse,dataLOO20)
dataLOO1AUC <- bootstrappingGLM(dataLOO1in,bootsamples)
FX20AUC<-mean(dataLOO1AUC)
quantiles20 <- quantile(dataLOO1AUC,probs=seq(0,1,.025), type=3,na.rm =TRUE)
print(FX20AUC)
print(quantiles20[quantilspostitions])
#varImpPlot(BED20,sort=TRUE,n.var=10,type=2)


dataLOO2in <- cbind(responses,binaryRepsonse,dataLOO40)
dataLOO2AUC <- bootstrappingGLM(dataLOO2in,bootsamples)
FX40AUC<-mean(dataLOO2AUC)
quantiles40 <- quantile(dataLOO2AUC,probs=seq(0,1,.025), type=3,na.rm =TRUE)
print(FX40AUC)
print(quantiles40[quantilspostitions])
#varImpPlot(BED40,sort=TRUE,n.var=10,type=2)


pdf("LiverBED20.pdf")
bed20mdlroc <- randomForest(x=dataLOO20,y=ptResponsesliver,ntree=500,replace=TRUE,importance=TRUE,mtry=2,header=TRUE)
bed40mdlroc<- randomForest(x=dataLOO40,y=ptResponsesliver,ntree=500,replace=TRUE,importance=TRUE,mtry=2,header=TRUE)
preds20 <-as.numeric(predict(bed20mdlroc,dataLOO20, type = "response"))
preds40 <-as.numeric(predict(bed40mdlroc,dataLOO40, type = "response"))
aucVector20<-NULL
aucVector40<-NULL
rocobj<-roc(response=Liverdata[,4],predictor=bed20mdlroc$votes[,2])

plot.roc(rocobj)
dataIN<-NULL
dataIN<-data.frame(responsed=Liverdata[,3],binaryresponse=Liverdata[,4],feature1=dataLOO20[,1],feature2=dataLOO20[,2])
check<-NULL
count<-0
pts<-1:length(ptResponsesliver)
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

pdf("LiverBED40.pdf")
rocobj40<-roc(response=Liverdata[,4],predictor=bed40mdlroc$votes[,2])
plot(rocobj40)
dataIN<-NULL
dataIN<-data.frame(responsed=Liverdata[,3],binaryresponse=Liverdata[,4],feature1=dataLOO40[,1],feature2=dataLOO40[,2])

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







#pancfeats <-c(6,37)
print(cor(fx1[c(6)],fx1[c(1)],method="spearman"))
print(cor(fx1[c(6)],fx1[c(18)],method="spearman"))
print(cor(fx1[c(37)],fx1[c(1)],method="spearman"))
print(cor(fx1[c(37)],fx1[c(23)],method="spearman"))




#print(cor(dataLOO[,1],dataLOO[,2],method="spearman"))
#Correlation test



# nmbrPts<-length(responses)
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
#   ptresponses <- Liverdata[ptstesting,4]
#   #print(ptresponses)
#   mdlDefault<-NULL
#   mdlDefault <- randomForest(x=dataLOO[ptstesting,],y=ptresponses,ntree=500,mtry=2)
#   predf <- predict(mdlDefault,dataLOO[pt_to_exclude,],type='response')
#   #print(predf)
#   pred[ii]<- predf
#   #print(pred)
#   respreal[ii]<- responses[pt_to_exclude]
#   #print(respreal)
# }
# rocobj<-roc(response=respreal,predictor=pred,ci="true")
# print(rocobj)
# plot(rocobj)
