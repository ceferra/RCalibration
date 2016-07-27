
setwd("/home/cesar/Dropbox/research/Jupyter/Rcal")

source("calfunc.R")
source("calibrationMethods.R")

set.seed(0) 
require(caret)

library(MASS)


#predprobcal<-PAV1vsAll(as.double(as.factor(as.double(trainY))),trainProbs,testProbs)
library("RWeka")


datos <- read.arff("data/kr-vs-kp.arff")
datos <- read.arff("data/credit-g.arff")

posParamEstudio<-length(datos[1,])
nomParamEstudio<-names(datos)[posParamEstudio]
#datos<-datos[complete.cases(datos),]
# nomclasse<-names(datos)[length(names(datos))]
tam<-length(datos[,1])
# print(tam)
#mnom<-nomclasses[id]
numclass<-(length(levels(datos[,posParamEstudio])))
ob<-paste(nomParamEstudio,"~.",sep="")
w<-datos
indiceStrat<-createDataPartition(w[,nomParamEstudio], p = 0.5,list=FALSE)
train<-w[indiceStrat,]
test<-w[-indiceStrat,] 

indiceStrat<-createDataPartition(train[,nomParamEstudio], p = 0.5,list=FALSE)
validation<-train[indiceStrat,]
trainv<-train[-indiceStrat,] 


 model <- J48(ob, data = trainv, control = Weka_control(U = TRUE,A=TRUE))
predprobtemp<-predict(model, newdata = test,type="prob")
predprobval<-predict(model, newdata = validation,type="prob")
#calibration
#predprobcal<-PAV1vsAll(as.factor(as.double(validation[[nomParamEstudio]])),predprobval,predprobtemp)
predprobcal<-plattCal(as.factor(as.double(validation[[nomParamEstudio]])),predprobval,predprobtemp)


NB <- make_Weka_classifier("weka/classifiers/bayes/NaiveBayes")
model <- NB(ob, data = trainv)
predprobtempNB<-predict(model, newdata = test,type="prob")
predprobvalNB<-predict(model, newdata = validation,type="prob")
#calibration
predprobcalNB<-plattCal(as.factor(as.double(validation[[nomParamEstudio]])),predprobvalNB,predprobtempNB)

testProbs <- data.frame(obs =test[[nomParamEstudio]] ,NB=predprobtempNB[,1],J48=predprobtemp[,1])

calPlotData <- calibration(obs ~ J48+NB, data = testProbs)

xyplot(calPlotData, auto.key = list(columns = 2))


res<-MSEdecomp(cbind(as.double(test[[nomParamEstudio]])-1,predprobtemp[,2]))
print(paste("J48 MSE",MSE2(cbind(as.double(test[[nomParamEstudio]])-1,predprobtemp[,2]))))
print(paste("J48 cal",res$cal,"ref",res$ref))

res<-MSEdecomp(cbind(as.double(test[[nomParamEstudio]])-1,predprobtempNB[,2]))
print(paste("NB MSE",MSE2(cbind(as.double(test[[nomParamEstudio]])-1,predprobtempNB[,2]))))
print(paste("NB cal",res$cal,"ref",res$ref))


testProbs <- data.frame(obs =test[[nomParamEstudio]] ,NB=predprobcalNB[,1],J48=predprobcal[,1])

calPlotData <- calibration(obs ~  J48+NB, data = testProbs)

xyplot(calPlotData, auto.key = list(columns = 2))

res<-MSEdecomp(cbind(as.double(test[[nomParamEstudio]])-1,predprobcal[,2]))
print(paste("calJ48 MSE",MSE2(cbind(as.double(test[[nomParamEstudio]])-1,predprobcal[,2]))))
print(paste("calJ48 cal",res$cal,"ref",res$ref))

res<-MSEdecomp(cbind(as.double(test[[nomParamEstudio]])-1,predprobcalNB[,2]))
print(paste("calNB MSE",MSE2(cbind(as.double(test[[nomParamEstudio]])-1,predprobcalNB[,2]))))
print(paste("calNB cal",res$cal,"ref",res$ref))
