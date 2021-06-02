# This program is written to analyze pancreas delta radiomics texture
# features as a function of delivered BED as follows:
#           Change in features should be calculated for change after delivery of
#           at least 20 Gy. EX: BED 20 Gy will be a delivered BED of 20 >= BED < 40
#                               BED 40 Gy will be a delivered BED of 40>= BED < 60
#                               BED 60 Gy will be a delivered BED of BED >= 60+ BED
# Written by Garrett Simpson gns24@miami.edu v1 3.9.2021
# Clear environment variables
# Command to copy to machine clipboard:
#       write.table(dataLOO, "clipboard", sep="\t", row.names=FALSE)

rm(list=ls()) 
# Call libraries used
library('randomForest')
library('pROC')
library('psych')
# Location of CSV with data organized as described above
PancData <- read.csv("C:\\Users\\gns24\\Documents\\PancE64_28pts_BED.csv",header=TRUE)
# Number of patients
npatients <- 1:nrow(data)
# Data is already organized by BED delivered per fraction for calculation of
# Texture feature changes as a function of dose
fx1 <- PancData[,5:43]
fx2 <- PancData[,44:82]
fx3 <- PancData[,83:121]
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
dataLOO <- BED20[c(BED20TopSix[r])]
#dataLOO <- BED40[c(BED40TopSix[r])]
responses <- PancData[,3]
nmbrPts<-dim(PancData)
pts<-1:nmbrPts[1]
pred<-NULL
responsepred <-NULL
for (i in pts){
  ptstesting<- pts[pts!=i]
  ptresponses <- PancData[pts!=i,3]
  mdlDefault <- randomForest(x=dataLOO[ptstesting,],y=ptresponses,ntree=500,mtry=2)
  pred[i] <- predict(mdlDefault,dataLOO[i,],type='response')
}
# trueResponse should be the two categories 'NR' or 'RS'
trueResponse <-responses
pred <- pred
rocobj<-roc(response=(trueResponse),predictor=pred,ci="true")
print((rocobj))
plot(rocobj)
#check for correlation between texture features
print(cor(dataLOO[,1],dataLOO[,2],method="spearman"))


