# This program is to test the different methods of deltaRadiomics 
# for pancreas paper fraction 1 features will be considered as calculated,
# features from fractions 2-5 will be normalized by their fraction 1 value.
# Two methods are under scrutiny. The first (1) is selecting 
# predictive features from each fraction using raw feature values and features
# normalized by fraction 1 values. The top features will also be determined for
# the delta features using normalized texture values and raw change in values.


# Data structure prior to running this code is of UTMOST importance. And
# should be structured as follows:
# 1st column: should patient IDs, 2nd: binary response in characters, 3rd: binary
# response in terms of 0/1, columns 4-42: texture features from fraction 1-> work by Garrett Simpson
# uses GLCM, GLRLM, GLSZM, and NGTDM features. As this work is to analyze SBRT features,
# there are only 5 sets of 39 columns of texture features. Explicitly, columns 43-81: fraction 2
# based features that are stored in the EXACT order of columns 4-42, Columns 82-120:
# is filled with fraction 3 features, columns 121-159: fraction 4, and columns 160-198:
# is fraction 5 features.


# clear environment variables
rm(list=ls()) 
# call libraries used
library('randomForest')
library('pROC')
library('psych')
library('boot')
library('caTools')
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

# location of CSV with data organized as described above
#data <- read.csv("C:\\Users\\gns24\\Documents\\PancE64_30pts_4.29.2021.csv",header=TRUE)
datap <- read.csv("C:\\Users\\gns24\\Documents\\PancE64_30pts_4.29.2021.csv",header=TRUE)

# number of patients
npatients <- 1:nrow(datap)
# organizing data for manipulation
# fx1data is the raw value of texture features from the csv file
fx1data <- datap[,4:42]
# fx2data is normalized by fraction 1 values
fx2data <- datap[,43:81]/fx1data
# ogfx2data is the value direct from the csv
ogfx2data<- datap[,43:81]
fx3data <- datap[,82:120]/fx1data
ogfx3data <- datap[,82:120]
fx4data <- datap[,121:159]/fx1data
ogfx4data <- datap[,121:159]
fx5data <- datap[,160:198]/fx1data
ogfx5data <- datap[,160:198]
# defining response variables and coding for use in randomforest algorithm
labels <- c("NR","RS")
levels <- c(0,1)
ptResponses <- factor(datap[,3],levels,labels)

  
# Define table variables for Gini importance scores for data analysis
a<-c("Features","Gini Score")
temp<-colnames(datap)
b<-temp[4:42]
txfx1names<-temp[4:42]
txfx2names<-temp[43:81]
txfx3names<-temp[82:120]
txfx4names<-temp[121:159]
txfx5names<-temp[160:198]


##Lines 63-73 no longer used. Can be converted to predict based upon fx1 values though
## Verify/check to see if average features hold up in the larger library.
## Original work selected GLSZM GLV and GLCM Energy
summedVals <- fx1data + ogfx2data + ogfx3data + ogfx4data + ogfx5data
avefeatures <- summedVals/5
aveFeatureMdl <- randomForest(x=avefeatures,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
aveimp <- importance(aveFeatureMdl)
aveorderImp <- order(aveimp[,4],decreasing=TRUE)
avetopsix <- txfx1names[aveorderImp[1:6]]
print(aveFeatureMdl)
print(avetopsix)
#bootstrapping performance
dataIN<-data.frame(responsed=ptResponses,binaryresponse=datap[,3],feature1=avefeatures[,3],feature2=avefeatures[,1])
nmbrPts<-length(dataIN[,1])
pts<-1:nmbrPts[1]
aucVector<-NULL
rocobj<-roc(response=datap[,2],predictor=aveFeatureMdl$votes[,2])
plot(rocobj)
for(bsn in 1:1000){
  sample<-NULL
  sample <- sample.split(pts, SplitRatio = 2/3)
  trainingpts<-NULL
  trainingpts<- subset(pts, sample == TRUE)
  LRModel<-NULL
  preds<-NULL
  LRModel <-glm(formula = responsed~feature1+feature2,data=dataIN[trainingpts,],family=binomial(logit), maxit = 100)
  preds <- predict.glm(LRModel,dataIN[,3:4],type="response")
  tmpPerform<-roc(dataIN[,1],preds)
  aucVector[bsn]<-as.numeric(auc(tmpPerform))
  plot.roc(tmpPerform,add=TRUE,col = '#D3D3D3')
}




rte<-c(33,37)
bootsamples<-1000
quantpositions <-c(2,40)
dataLOOAver <- avefeatures[rte]
dataLOOavein <- cbind(responses,datap[,3],dataLOOAver)
dataLOOaveAUC <- bootstrappingGLM(dataLOOavein,bootsamples)
FXAVEAUC<-mean(dataLOOaveAUC)
quantilesave <- quantile(dataLOOaveAUC,probs=seq(0,1,.025), type=3,na.rm =TRUE)
print(FXAVEAUC)
print(quantilesave[quantpositions])
varImpPlot(aveFeatureMdl,sort=TRUE,n.var=10,type=2,main = c("Variable Importance Plot"))


# Normalized feature importance estimations
# Models trained and importances calculated
fx1FeatureMdl <- randomForest(x=fx1data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
fx1imp <- importance(fx1FeatureMdl)
fx1orderImp <- order(fx1imp[,4],decreasing=TRUE)
fx1topsix <- txfx1names[fx1orderImp[1:6]]
fx2FeatureMdl <- randomForest(x=fx2data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
fx2imp <- importance(fx2FeatureMdl)
fx2orderImp <- order(fx2imp[,4],decreasing=TRUE)
fx2topsix <- txfx2names[fx2orderImp[1:6]]
fx3FeatureMdl <- randomForest(x=fx3data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
fx3imp <- importance(fx3FeatureMdl)
fx3orderImp <- order(fx3imp[,4],decreasing=TRUE)
fx3topsix <- txfx3names[fx3orderImp[1:6]]
fx4FeatureMdl <- randomForest(x=fx4data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
fx4imp <- importance(fx4FeatureMdl)
fx4orderImp <- order(fx4imp[,4],decreasing=TRUE)
fx4topsix <- txfx4names[fx4orderImp[1:6]]
fx5FeatureMdl <- randomForest(x=fx5data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
fx5imp <- importance(fx5FeatureMdl)
fx5orderImp <- order(fx5imp[,4],decreasing=TRUE)
fx5topsix <- txfx5names[fx5orderImp[1:6]]

# Most important features in one place
# Extract top 6 features (ie, mtry value, sqrt(39))

selectedFeaturesNorm<-cbind(fx1topsix,fx2topsix)
selectedFeaturesNorm<-cbind(selectedFeaturesNorm,fx3topsix)
selectedFeaturesNorm<-cbind(selectedFeaturesNorm,fx4topsix)
selectedFeaturesNorm<-cbind(selectedFeaturesNorm,fx5topsix)

#write.csv(selectedFeaturesNorm,"C:\\Users\\gns24\\Desktop\\TA Studies\\Panc DeltaTA update\\DeltaRadTesting25pts\\pancreasRFABSw28features.csv", row.names = FALSE)


# Normalized delta radiomics features for each fraction
# produce normalized delta models for fx 1-2/1, fx 1-3/1, fx 1-4/1, and fx 1-5/1
# Predefining some of the variables

temp<-colnames(datap)
b<-temp[4:198]

fx1data <- datap[,4:42]
ogfx2data<- datap[,43:81]
ogfx3data <- datap[,82:120]
ogfx4data <- datap[,121:159]
ogfx5data <- datap[,160:198]

# Normalized delta radiomics texture features
fx12data<-(fx1data-ogfx2data)/fx1data
fx13data<-(fx1data-ogfx3data)/fx1data
fx14data<-(fx1data-ogfx4data)/fx1data
fx15data<-(fx1data-ogfx5data)/fx1data

# Train normalized delta radiomics models and find important values for
# all of the fractional changes
fx2deltaNorm <- randomForest(x=fx12data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
fx2deltaNormimp <- importance(fx2deltaNorm)
fx2deltaNormorderImp <- order(fx2deltaNormimp[,4],decreasing=TRUE)
fx2cumulativetopsix <- b[fx2deltaNormorderImp[1:6]]
fx3deltaNorm <- randomForest(x=fx13data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
fx3deltaNormimp <- importance(fx3deltaNorm)
fx3deltaNormimpImp <- order(fx3deltaNormimp[,4],decreasing=TRUE)
fx3cumulativetopsix <- b[fx3deltaNormimpImp[1:6]]
fx4deltaNorm <- randomForest(x=fx14data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
fx4deltaNormimp <- importance(fx4deltaNorm)
fx4deltaNormorderImp <- order(fx4deltaNormimp[,4],decreasing=TRUE)
fx4cumulativetopsix <- b[fx4deltaNormorderImp[1:6]]
fx5deltaNorm <- randomForest(x=fx15data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
fx5deltaNormimp <- importance(fx5deltaNorm)
fx5deltaNormorderImp <- order(fx5deltaNormimp[,4],decreasing=TRUE)
fx5cumulativetopsix <- b[fx5deltaNormorderImp[1:6]]

fx2deltaNormtopsixp<- c(fx2cumulativetopsix)
fx3deltaNormtopsixp<- c(fx3cumulativetopsix)
fx4deltaNormtopsixp<- c(fx4cumulativetopsix)
fx5deltaNormtopsixp<- c(fx5cumulativetopsix)

selectedFeaturesDelta<-cbind(fx2deltaNormtopsixp,fx3deltaNormtopsixp)
selectedFeaturesDelta<-cbind(selectedFeaturesDelta,fx4deltaNormtopsixp)
selectedFeaturesDelta<-cbind(selectedFeaturesDelta,fx5deltaNormtopsixp)


# Output information to csvfile for optimal viewing
# Features being selected should be used for delta radiomics analysis
# Change file name based upon the type of features used
write.csv(selectedFeaturesDelta,"C:\\Users\\gns24\\Desktop\\TA Studies\\Panc DeltaTA update\\DeltaRadTesting25pts\\pancreasRFNormDeltaFeatures28Delta.csv", row.names = FALSE)

# Selection of features based upon raw values (no normalization by fx1 values)
# delta radiomics texture features
Rawfx12data<-(fx1data-ogfx2data)
Rawfx13data<-(fx1data-ogfx3data)
Rawfx14data<-(fx1data-ogfx4data)
Rawfx15data<-(fx1data-ogfx5data)

#Rawfx12data<-(fx1data-ogfx2data)/fx1data
#Rawfx13data<-(fx1data-ogfx3data)/fx1data
#Rawfx14data<-(fx1data-ogfx4data)/fx1data
#Rawfx15data<-(fx1data-ogfx5data)/fx1data


# Train models and rank importance for absolute change between fractions for change
# between fraction 1 and fractions 1-4
Rawfx2delta <- randomForest(x=Rawfx12data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
Rawfx2deltaimp <- importance(Rawfx2delta)
Rawfx2deltaorderImp <- order(Rawfx2deltaimp[,4],decreasing=TRUE)
Rawfx2cumulativetopsix <- b[Rawfx2deltaorderImp[1:6]]
#varImpPlot(Rawfx2delta,sort=TRUE,n.var=10,type=2)

Rawfx3delta <- randomForest(x=Rawfx13data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=61,header=TRUE)
Rawfx3deltaimp <- importance(Rawfx3delta)
Rawfx3deltaimpImp <- order(Rawfx3deltaimp[,4],decreasing=TRUE)
Rawfx3cumulativetopsix <- b[Rawfx3deltaimpImp[1:6]]
#varImpPlot(Rawfx3delta,sort=TRUE,n.var=10,type=2)

Rawfx4delta <- randomForest(x=Rawfx14data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
Rawfx4deltaimp <- importance(Rawfx4delta)
Rawfx4deltaorderImp <- order(Rawfx4deltaimp[,4],decreasing=TRUE)
Rawfx4cumulativetopsix <- b[Rawfx4deltaorderImp[1:6]]
#varImpPlot(Rawfx4delta,sort=TRUE,n.var=10,type=2)

Rawfx5delta <- randomForest(x=Rawfx15data,y=ptResponses,ntree=500,replace=TRUE,importance=TRUE,mtry=6,header=TRUE)
Rawfx5deltaimp <- importance(Rawfx5delta)
Rawfx5deltaorderImp <- order(Rawfx5deltaimp[,4],decreasing=TRUE)
Rawfx5cumulativetopsix <- b[Rawfx5deltaorderImp[1:6]]
#varImpPlot(Rawfx5delta,sort=TRUE,n.var=10,type=2)

Rawfx2deltatopsixp<- c(Rawfx2cumulativetopsix)
Rawfx3deltatopsixp<- c(Rawfx3cumulativetopsix)
Rawfx4deltatopsixp<- c(Rawfx4cumulativetopsix)
Rawfx5deltatopsixp<- c(Rawfx5cumulativetopsix)

RawselectedFeaturesDelta<-cbind(Rawfx2deltatopsixp,Rawfx3deltatopsixp)
RawselectedFeaturesDelta<-cbind(RawselectedFeaturesDelta,Rawfx4deltatopsixp)
RawselectedFeaturesDelta<-cbind(RawselectedFeaturesDelta,Rawfx5deltatopsixp)

print(RawselectedFeaturesDelta)

#write.csv(RawselectedFeaturesDelta,"C:\\Users\\gns24\\Desktop\\TA Studies\\Panc DeltaTA update\\DeltaRadTesting25pts\\pancreasRFAbsValDeltaFeatures28Delta.csv", row.names = FALSE)


# This section is to test the predictive modeling using leave-n-out, n=1,
# will use values selected from previous section. 
# fx12data is normalized by the first fraction value
# Rawfx12data is raw value w/o normalization
# avefeatures is the averaged features

#selecting which columns (features to use)
r<-c(1,2)
#Select which type of features you want to run using the leave-one-out
#reference top 2 features automatically by variable name
#Normalized features are in fx12data and absolute change values are in Rawfx12data
#dataLOO <- fx12data[c(fx2cumulativetopsix[r])]
dataLOO <- Rawfx12data[c(Rawfx2cumulativetopsix[r])]
#dataLOO <- fx12data[,r]
#dataLOO <- avefeatures[c(avetopsix[r])]
responses <- datap[,2]
nmbrPts<-dim(avefeatures)
pts<-1:nmbrPts[1]
pred<-NULL
responsepred <-NULL
#for (i in pts){
  #ptstesting<- pts[pts!=i]
  #ptresponses <- data[pts!=i,3]
  #mdlDefault <- randomForest(x=dataLOO[ptstesting,],y=ptresponses,ntree=500,mtry=2)
  #pred[i] <- predict(mdlDefault,dataLOO[i,],type='response')
#}
bootsamples<-1000
quantpositions <-c(2,40)
dataLOO1 <- Rawfx12data[c(Rawfx2cumulativetopsix[r])]
dataLOO1in <- cbind(responses,datap[,3],dataLOO1)
dataLOO1AUC <- bootstrappingGLM(dataLOO1in,bootsamples)
FX12AUC<-mean(dataLOO1AUC)
quantiles12 <- quantile(dataLOO1AUC,probs=seq(0,1,.025), type=3,na.rm =TRUE)
print(FX12AUC)
print(quantiles12[quantpositions])

dataLOO2 <- Rawfx13data[c(Rawfx3cumulativetopsix[r])]
dataLOO2in <- cbind(responses,datap[,3],dataLOO2)
dataLOO2AUC <- bootstrappingGLM(dataLOO2in,bootsamples)
FX13AUC<-mean(dataLOO2AUC)
quantiles13 <- quantile(dataLOO2AUC,probs=seq(0,1,.025), type=3,na.rm =TRUE)
print(FX13AUC)
print(quantiles13[quantpositions])

dataLOO3 <- Rawfx14data[c(Rawfx4cumulativetopsix[r])]
dataLOO3in <- cbind(responses,datap[,3],dataLOO3)
dataLOO3AUC <- bootstrappingGLM(dataLOO3in,bootsamples)
FX14AUC<-mean(dataLOO3AUC)
quantiles14 <- quantile(dataLOO3AUC,probs=seq(0,1,.025), type=3,na.rm =TRUE)
print(FX14AUC)
print(quantiles14[quantpositions])

dataLOO4 <- Rawfx15data[c(Rawfx5cumulativetopsix[r])]
dataLOO4in <- cbind(responses,datap[,3],dataLOO4)
dataLOO4AUC <- bootstrappingGLM(dataLOO4in,bootsamples)
FX15AUC<-mean(dataLOO4AUC)
quantiles15 <- quantile(dataLOO4AUC,probs=seq(0,1,.025), type=3,na.rm =TRUE)
print(FX15AUC)
print(quantiles15[quantpositions])


# trueResponse should be the two categories 'NR' or 'RS'
trueResponse <-responses
pred <- pred
rocobj<-roc(response=(trueResponse),predictor=pred,ci="true")
print((rocobj))
plot(rocobj)
#check for correlation between texture features
print(cor(dataLOO[,1],dataLOO[,2],method="spearman"))

print(cor(fx20BED[c(6)],fx20BED[c(BED20Topsix[1])],method="spearman"))
print(cor(fx20BED[c(6)],fx20BED[c(BED20Topsix[2])],method="spearman"))
print(cor(fx20BED[c(37)],fx20BED[c(BED20Topsix[1])],method="spearman"))
print(cor(fx20BED[c(37)],fx20BED[c(BED20Topsix[2])],method="spearman"))



#The following are two examples of different box plots to compare features between responders and non-responders
# plot(data[,2],fx12data[,31],xlab="Binary Patient Response",ylab="Absolute Change: GLSZM LZLGL")
# plot(data[,2],Rawfx12data[,6],xlab="Binary Patient Response",ylab="Absolute Change: GLCM Sum Average")


##The following is the base code to perform  hierarchical cluster analysis with the data
#library('cluster')
##data for clustering must be standardized (mean<-0,sd<-1)
#standTextureFeatures <-scale(cbind(data[,2],fx12data))
##divisive cluster
#hc1 <- agnes(standTextureFeatures,method='complete')
#pdf("clusterAnalysis.pdf")
#plot(hc1)
#dev.off()
#agnes computes the agglomerative coefficient, which measures the amount of clustering
## structure found (values closer to 1 suggest strong clustering structure). 
#hc1$ac


# The follwoing code is written to find any correlations between the most predictive features
# selected from initial analysis.
# This retieves texture names from the data instance used for training the RF
# 
correlationChecks<-fx12data[c(Rawfx2cumulativetopsix)]
#Option to save as pdf file
pdf("clusterAnalysis.pdf")
pairs.panels(correlationChecks)
dev.off()

#correlation between GRLRM and GLSZM

print(cor.test(fx1data[,c(9)],fx1data[,c(22)],alternative = c("two.sided"), conf.level = 0.95,method="spearman"))
print(cor.test(fx1data[,c(10)],fx1data[,c(23)],method="spearman"))
print(cor.test(fx1data[,c(11)],fx1data[,c(24)],method="spearman"))
print(cor.test(fx1data[,c(12)],fx1data[,c(25)],method="spearman"))
print(cor.test(fx1data[,c(13)],fx1data[,c(26)],method="spearman"))
print(cor.test(fx1data[,c(14)],fx1data[,c(27)],method="spearman"))
print(cor.test(fx1data[,c(15)],fx1data[,c(28)],method="spearman"))
print(cor.test(fx1data[,c(16)],fx1data[,c(29)],method="spearman"))
print(cor.test(fx1data[,c(17)],fx1data[,c(30)],method="spearman"))
print(cor.test(fx1data[,c(18)],fx1data[,c(31)],method="spearman"))
print(cor.test(fx1data[,c(19)],fx1data[,c(32)],method="spearman"))
print(cor.test(fx1data[,c(20)],fx1data[,c(33)],method="spearman"))
print(cor.test(fx1data[,c(21)],fx1data[,c(34)],method="spearman"))







