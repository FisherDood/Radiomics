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