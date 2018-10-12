set.seed(123)

#Aqueous solubility
aqsolubil.data<-data%>%filter(description=='Aqueous solubility')
aqsolubil.data<-aqsolubil.data%>%filter(published_units!='mg ml-1')#só nos nulls

#remove NA
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(mw_freebase))
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(psa))
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(acd_most_apka))
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(hba))
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(hbd))
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(acd_logp))
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(acd_logd))
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(rtb))
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(acd_most_bpka))
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(mw_monoisotopic))
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(aromatic_rings))
aqsolubil.data<-aqsolubil.data%>%filter(!is.na(full_mwt))

aqsolubil.data<-aqsolubil.data[!duplicated(aqsolubil.data),]

###training

samp <- sample(nrow(aqsolubil.data), 0.6 * nrow(aqsolubil.data))
aqsolubil.train <- aqsolubil.data[samp, ]#744
aqsolubil.test <- aqsolubil.data[-samp, ]#497

aqsolubil.x.train<-aqsolubil.train%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,rtb,mw_monoisotopic, aromatic_rings,full_mwt) # 68%
aqsolubil.y.train<-aqsolubil.train$standard_value

################################random forest################################
aqsolubil.rf.mdl <-randomForest(aqsolubil.y.train,x=aqsolubil.x.train,mtry = 7,ntree = 1500)
print(aqsolubil.rf.mdl)
varImpPlot(aqsolubil.rf.mdl)
importance(aqsolubil.rf.mdl)


aqsolubil.pred<-predict(aqsolubil.rf.mdl,newdata=aqsolubil.test)
plot(aqsolubil.pred,aqsolubil.test$standard_value)
abline(0,1,col="blue")

###MSE formulation
rf.PredictionRMSE<-mean((aqsolubil.pred-aqsolubil.test$standard_value)^2) #1.23 
rf.PredictionMAE<-mean(abs(aqsolubil.pred-aqsolubil.test$standard_value))

#tuning rf
aqsolubil.rf.mdl.tuned<-tuneRF(aqsolubil.x.train,aqsolubil.y.train,mtryStart = 3,ntreeTry = 1500, improve = 0.01,doBest = TRUE)
aqsolubil.pred.tuned<-predict(aqsolubil.rf.mdl.tuned,aqsolubil.test)
plot(aqsolubil.pred.tuned,aqsolubil.test$standard_value)
abline(0,1, col="blue")
rf.Prediction.tunedRMSE<-mean((aqsolubil.pred.tuned-aqsolubil.test$standard_value)^2) 
print(rf.Prediction.tunedRMSE)#1.19
rf.Prediction.tunedMAE<-mean(abs(aqsolubil.pred.tuned-aqsolubil.test$standard_value))
print(rf.Prediction.tunedMAE)#0.59

######################################SVR#######################################

aqsolubil.x.train<-aqsolubil.train%>%select(full_mwt,mw_freebase,psa,acd_most_apka,hba,hbd,acd_logp,acd_logd,rtb,acd_most_bpka,mw_monoisotopic,aromatic_rings) # 68%

aqsolubil.modelsvr<-svm(y=aqsolubil.y.train,x=aqsolubil.x.train)
aqsolubil.predict.svr <-predict(aqsolubil.modelsvr,aqsolubil.x.train)

print(aqsolubil.modelsvr)

plot(aqsolubil.train$standard_value)
points(aqsolubil.predict.svr,col="Green")
aqsolubil.modelsvr.error <- aqsolubil.train$standard_value-aqsolubil.predict.svr
aqsolubil.svrPredictionRMSE <-sqrt(mean(aqsolubil.modelsvr.error^2))
print (aqsolubil.svrPredictionRMSE)#1.14
aqsolubil.svrPredictionMAE<-mean(abs(aqsolubil.predict.svr-aqsolubil.train$standard_value))


tc <- tune.control(cross = 10)

aqsolubil.svr.tunedresult<-tune(svm,train.x=aqsolubil.x.train,train.y=aqsolubil.y.train,ranges=list(epsilon=seq(0,1,0.1),cost=2^(2:9)))
print(aqsolubil.svr.tunedresult)

#epsilon cost
#0.1   32 
aqsolubil.svr.tunedresult<-tune(svm,train.x=aqsolubil.x.train,train.y=aqsolubil.y.train,ranges=list(epsilon=seq(0,0.2,0.01),cost=2^(2:7)),tunecontrol = tc)

print(aqsolubil.svr.tunedresult)
plot(aqsolubil.svr.tunedresult)

aqsolubil.predict.svr.tuned <-predict(aqsolubil.svr.tunedresult$best.model,aqsolubil.x.train)
aqsolubil.modelsvr.error.tuned <- aqsolubil.train$standard_value-aqsolubil.predict.svr.tuned
svrPrediction.tunedRMSE <-sqrt(mean(aqsolubil.modelsvr.error.tuned^2))
print (svrPrediction.tunedRMSE)#0.62
svrPrediction.tunedMAE <- mean(abs(aqsolubil.predict.svr.tuned-aqsolubil.train$standard_value))
print(svrPrediction.tunedMAE)


################################BAYES############################################
samp <- sample(nrow(aqsolubil.data), 0.6 * nrow(aqsolubil.data))
aqsolubil.train <- aqsolubil.data[samp, ]#744
aqsolubil.test <- aqsolubil.data[-samp, ]#497

aqsolubil.train.net<-aqsolubil.train%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,rtb,mw_monoisotopic, aromatic_rings,full_mwt,standard_value) # 76%
aqsolubil.test.net<-aqsolubil.test%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,rtb,mw_monoisotopic, aromatic_rings,full_mwt,standard_value) # 76%


aqsolubil.hc.net <- hc(aqsolubil.train.net, debug=TRUE)#hill climbing

# plot network
plot(aqsolubil.hc.net)

# fit the network with the data itself
aqsolubil.hc.fit <- bn.fit(aqsolubil.hc.net, aqsolubil.train.net)

aqsolubil.hc.predict <-predict(aqsolubil.hc.fit,"standard_value",aqsolubil.test.net)

plot(aqsolubil.hc.predict)
points(aqsolubil.test.net$standard_value,col="Red")

#error
aqsolubil.modelhc.error <- aqsolubil.test.net$standard_value-aqsolubil.hc.predict
aqsolubil.hc.PredictionRMSE <-sqrt(mean(aqsolubil.modelhc.error^2))
print (aqsolubil.hc.PredictionRMSE)#1.91
aqsolubil.hc.PredictionMAE<-mean(abs(aqsolubil.hc.predict-aqsolubil.test.net$standard_value))

cbind(predicted=aqsolubil.hc.predict,actual=aqsolubil.test.net$standard_value)


###########################Multiple Linear regression##############################
plot(log(aqsolubil.train$standard_value+11))

aqsolubil.lin.reg <- lm(log(standard_value+11) ~ mw_freebase +  psa + acd_most_apka + acd_most_bpka + hba +
                          hbd + rtb+mw_monoisotopic+aromatic_rings+full_mwt, data = aqsolubil.train)


# Inspect the model
summary(aqsolubil.lin.reg)

aqsolubil.pred.lin <- exp(predict(aqsolubil.lin.reg,aqsolubil.test))-11

aqsolubil.lin.reg.RMSE <- sqrt(mean((aqsolubil.pred.lin-aqsolubil.test$standard_value)^2))
print(aqsolubil.lin.reg.RMSE)
aqsolubil.lin.reg.MAE <- mean(abs(aqsolubil.pred.lin-aqsolubil.test$standard_value))
print(aqsolubil.lin.reg.MAE)

################################# Neural Networks ###################################

aqsolubil.data.NN<-aqsolubil.data%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,rtb,mw_monoisotopic, aromatic_rings,full_mwt,standard_value) # 76%

#check if no na is available
apply(aqsolubil.data.NN,2,function(x) sum(is.na(x)))#fun?

aqsolubil.lm.fit<-glm(standard_value~.,data=aqsolubil.train.net)

summary(aqsolubil.lm.fit)

aqsolubil.pr.lm<-predict(aqsolubil.lm.fit,aqsolubil.test.net)
aqsolubil.MSE.lm<-sum((aqsolubil.pr.lm-aqsolubil.test.net$standard_value)^2)/nrow(aqsolubil.test)
print(aqsolubil.MSE.lm)



aqsolubil.maxs<-apply(aqsolubil.data.NN,2,max)
aqsolubil.mins<-apply(aqsolubil.data.NN,2,min)

aqsolubil.scaled<-as.data.frame(scale(aqsolubil.data.NN,center=aqsolubil.mins,scale=aqsolubil.maxs-aqsolubil.mins))

aqsolubil.train.NN<-aqsolubil.scaled[samp,]
aqsolubil.test.NN<-aqsolubil.scaled[-samp,]

aqsolubil.n<-names(aqsolubil.train.NN)
aqsolubil.f<-as.formula(paste("standard_value ~",paste(aqsolubil.n[!aqsolubil.n %in% "standard_value"],collapse = " + ")))

aqsolubil.NN.model<-neuralnet(aqsolubil.f,data=aqsolubil.train.NN,hidden=c(6,3),linear.output = TRUE)
#algorithm did not converge in 1 of 1 repetition(s) within the stepmax 

plot(aqsolubil.NN.model)

str(aqsolubil.test.NN[,1:10])#retirar o standard_value
aqsolubil.pred.NN<-compute(aqsolubil.NN.model,aqsolubil.test.NN[,1:10])

aqsolubil.pred.NN_<-aqsolubil.pred.NN$net.result*(max(aqsolubil.data.NN$standard_value)-min(aqsolubil.data.NN$standard_value))+min(aqsolubil.data.NN$standard_value)
aqsolubil.pred.NN.response<-(aqsolubil.test.NN$standard_value)*(max(aqsolubil.data.NN$standard_value)-min(aqsolubil.data.NN$standard_value))+min(aqsolubil.data.NN$standard_value)

aqsolubil.nn.RMSE<-sum((aqsolubil.pred.NN.response-aqsolubil.pred.NN_)^2)/nrow(aqsolubil.test.NN)
print(aqsolubil.nn.RMSE)

aqsolubil.nn.MAE<-mean(abs(aqsolubil.pred.NN.response-aqsolubil.pred.NN_))
print(aqsolubil.nn.MAE)


plot(aqsolubil.test.NN$standard_value,aqsolubil.pred.NN$net.result,col="blue")
abline(0,1,lwd=2)


#function for cross-validation Neural Network.
aqsolubil.cv.NN.error<-NULL
k<-10

pbar <- create_progress_bar('text')
pbar$init(k)

for(i in 1:k){
  index <- sample(1:nrow(aqsolubil.data),round(0.9*nrow(aqsolubil.data)))
  aqsolubil.train.NN.cv <- aqsolubil.scaled[index,]
  aqsolubil.test.NN.cv <- aqsolubil.scaled[-index,]
  
  aqsolubil.nn.cv <- neuralnet(aqsolubil.f,data=aqsolubil.train.NN.cv,hidden=c(6,2),linear.output=T)
  
  aqsolubil.pr.nn.cv <- compute(aqsolubil.nn.cv,aqsolubil.test.NN.cv[,1:10])
  aqsolubil.pr.nn.cv <- aqsolubil.pr.nn.cv$net.result*(max(aqsolubil.data$standard_value)-min(aqsolubil.data$standard_value))+min(aqsolubil.data$standard_value)
  
  aqsolubil.test.cv.r <- (aqsolubil.test.NN.cv$standard_value)*(max(aqsolubil.data$standard_value)-min(aqsolubil.data$standard_value))+min(aqsolubil.data$standard_value)
  
  aqsolubil.cv.NN.error[i] <- sum((aqsolubil.test.cv.r - aqsolubil.pr.nn.cv)^2)/nrow(aqsolubil.test.NN.cv)
  
  pbar$step()
}

mean(aqsolubil.cv.NN.error)

boxplot(aqsolubil.cv.NN.error,xlab='MSE NN CV',col='green',
        border='blue',names='CV error (MSE)',
        main='CV error (MSE) for NN',horizontal=TRUE)


######################################TOTAL#######################################
aqsolubil.accuracy <- data.frame(Method = c("Random Forest","SVR","HC.Bayes","SVR Tuned","Random forest Tuned","Multiple Linear Regression","Neural Networks"),
                       RMSE   = c(rf.PredictionRMSE,aqsolubil.svrPredictionRMSE,aqsolubil.hc.PredictionRMSE,svrPrediction.tunedRMSE,rf.Prediction.tunedRMSE,aqsolubil.lin.reg.RMSE,aqsolubil.nn.RMSE),
                       MAE    = c(rf.PredictionMAE,aqsolubil.svrPredictionMAE,aqsolubil.hc.PredictionMAE,svrPrediction.tunedMAE,rf.Prediction.tunedMAE,aqsolubil.lin.reg.MAE,aqsolubil.nn.MAE)) 


# Round the values and print the table
aqsolubil.accuracy$RMSE <- round(aqsolubil.accuracy$RMSE,2)
aqsolubil.accuracy$MAE <- round(aqsolubil.accuracy$MAE,2) 

print (aqsolubil.accuracy)
grid.table(aqsolubil.accuracy)
png("images/aqsolubil.png", height = 40*nrow(aqsolubil.accuracy), width = 100*ncol(aqsolubil.accuracy))
grid.table(aqsolubil.accuracy)
dev.off()






######################################OTHER#######################################WIP




###rf Cross validation 2
aqsolubil.rf.cv <- rf.crossValidation(aqsolubil.rf.mdl, aqsolubil.x.train, p=0.60, n=10, ntree=1500)# segunda versao cv
print(aqsolubil.rf.cv)

aqsolubil.rf.cv.tuned <- rf.crossValidation(aqsolubil.rf.mdl.tuned, aqsolubil.x.train, p=0.10, n=10, ntree=1500)# segunda versao cv
print(aqsolubil.rf.cv.tuned)



############CV TOTAL ##########
aqsolubil.accuracy.cv<-data.frame(Method=c("Random Forest CV10","Random Forest Tuned CV10","Neural Networks CV20","SVR Tuned CV10")
                             ,RMSE=c(mean(aqsolubil.rf.cv$y.rmse),mean(aqsolubil.rf.cv.tuned$y.rmse),mean(aqsolubil.cv.NN.error),svrPrediction.tunedRMSE)
                             ,MAE=c(mean(aqsolubil.rf.cv$y.mae),mean(aqsolubil.rf.cv.tuned$y.mae),0,svrPrediction.tunedMAE))

aqsolubil.accuracy.cv$RMSE <- round(aqsolubil.accuracy.cv$RMSE,2)
aqsolubil.accuracy.cv$MAE <- round(aqsolubil.accuracy.cv$MAE,2) 


print(aqsolubil.accuracy.cv)


