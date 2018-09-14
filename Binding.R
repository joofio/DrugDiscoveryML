set.seed(123)

#select assay necessary
binding.data<-data%>%filter(assay_id.1==225074)#Binding energy by using the equation deltaG obsd = -RT ln KD
print(assaydata)
summary(assaydata)

#eliminate NAs
binding.data<-binding.data%>%filter(!is.na(mw_freebase))
binding.data<-binding.data%>%filter(!is.na(psa))
binding.data<-binding.data%>%filter(!is.na(acd_most_apka))#quase 2000
binding.data<-binding.data%>%filter(!is.na(hba))
binding.data<-binding.data%>%filter(!is.na(hbd))
binding.data<-binding.data%>%filter(!is.na(acd_logp))
binding.data<-binding.data%>%filter(!is.na(acd_logd))
binding.data<-binding.data%>%filter(!is.na(rtb))
binding.data<-binding.data%>%filter(!is.na(acd_most_bpka))#quase 1000
binding.data<-binding.data%>%filter(!is.na(mw_monoisotopic))
binding.data<-binding.data%>%filter(!is.na(aromatic_rings))
binding.data<-binding.data%>%filter(!is.na(full_mwt))
binding.data<-binding.data%>%filter(!is.na(qed_weighted))


#tentar todas as possiveis para avaliar impt
binding.response<-binding.data$standard_value




binding.predictors1<-binding.data%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,acd_logp,acd_logd,heavy_atoms,rtb,hba_lipinski,mw_monoisotopic, molecular_species, aromatic_rings,full_mwt) # 64%
binding.predictors2<-binding.data%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,acd_logp,acd_logd,heavy_atoms,rtb,hba_lipinski,mw_monoisotopic, molecular_species, aromatic_rings,full_mwt,mw_monoisotopic,hbd_lipinski) # 64%
binding.predictors3<-binding.data%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,acd_logp,acd_logd,heavy_atoms) # 57%
binding.predictors4<-binding.data%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,acd_logp,acd_logd,qed_weighted) # 59%
binding.predictors5<-binding.data%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,acd_logp,acd_logd,qed_weighted) # 61%

################################random forest################################
###para Cross validation
binding.rf.mdl <-randomForest(binding.response,x=binding.predictors2)
print(binding.rf.mdl)
varImpPlot(binding.rf.mdl)
importance(binding.rf.mdl)
getTree(binding.rf.mdl, 3, labelVar=TRUE)#print de arvores


###para test e train
samp <- sample(nrow(binding.data), 0.6 * nrow(binding.data))
binding.train <- binding.data[samp, ]
binding.test <- binding.data[-samp, ]

binding.y.train<-binding.train$standard_value
binding.x.train<-binding.train%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,acd_logp,acd_logd,qed_weighted) # 59%
binding.x.train2<-binding.train%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,acd_logp,acd_logd,heavy_atoms,rtb,hba_lipinski,mw_monoisotopic, molecular_species, aromatic_rings,full_mwt,mw_monoisotopic,hbd_lipinski) # 64%

binding.rf.mdl2 <-randomForest(binding.y.train,x=binding.x.train2,ntree = 1500,proximity = TRUE )
print(binding.rf.mdl2)
binding.pred<-predict(binding.rf.mdl2,binding.test)


#MSE
plot(binding.pred,binding.test$standard_value)
abline(0,1, col="blue")

###MSE propriamente dito, e n previsto como no print rf
rf.predction.RMSE<-mean((binding.pred-binding.test$standard_value)^2)  #3.19
print(rf.predction.RMSE)
rf.predction.MAE<-mean(abs(binding.pred-binding.test$standard_value))
print(rf.predction.MAE)

#########?????MDSplot####TUNERF#####
help("MDSplot")
MDSplot(binding.rf.mdl2,fac=binding.y.train,pallete=c(1,15),pch=20)


#tuning rf
binding.rf.mdl.tuned<-tuneRF(binding.x.train2,binding.y.train,mtryStart = 3,ntreeTry = 1500, improve = 0.01,doBest = TRUE)
binding.pred.tuned<-predict(binding.rf.mdl.tuned,binding.test)
plot(binding.pred.tuned,binding.test$standard_value)
abline(0,1, col="blue")
binding.rf.Prediction.tunedRMSE<-mean((binding.pred.tuned-binding.test$standard_value)^2)  #2.77
print(binding.rf.Prediction.tunedRMSE)
binding.rf.Prediction.tunedMAE<-mean(abs(binding.pred.tuned-binding.test$standard_value))
print(binding.rf.Prediction.tunedMAE)#1.04

print(binding.rf.mdl.tuned)
######################################SVR#######################################

binding.modelsvr<-svm(y=binding.y.train,x=binding.x.train)
binding.predict.svr <-predict(binding.modelsvr,binding.x.train)

print(binding.modelsvr)

plot(binding.train$standard_value,binding.predict.svr)

binding.modelsvr.error <- binding.train$standard_value-binding.predict.svr
binding.svrPredictionRMSE <-sqrt(mean(binding.modelsvr.error^2))
print (binding.svrPredictionRMSE)#1.93
binding.svrPredictionMAE<-mean(abs(binding.predict.svr-binding.train$standard_value))
print(binding.svrPredictionMAE)#1.91

#tuning
binding.svr.tunedresult<-tune(svm,train.x=binding.x.train,train.y=binding.y.train,ranges=list(epsilon=seq(0,1,0.1),cost=2^(2:9)))#1.25
print(binding.svr.tunedresult)

#epsilon cost
#0.1   32 
#voltar a correr o tunedrsult com o seq de 0-0.1, com 0.01 e cost aumentar
tc <- tune.control(cross = 10)

binding.svr.tunedresult<-tune(svm,train.x=binding.x.train,train.y=binding.y.train,ranges=list(epsilon=seq(0,0.1,0.01),cost=2^(2:7)),tunecontrol=tc)#1.43

binding.svr.tunedresult<-tune(svm,train.x=binding.x.train,train.y=binding.y.train,ranges=list(epsilon=seq(0,0.4,0.025),cost=2^(2:7)),tunecontrol=tc)#1.22


print(binding.svr.tunedresult)
plot(binding.svr.tunedresult)

binding.predict.svr.tuned <-predict(binding.svr.tunedresult$best.model,binding.x.train)
binding.predict.svr.tuned <-predict(binding.svr.tunedresult$best.model,binding.test)

binding.modelsvr.error.tuned <- binding.train$standard_value-binding.predict.svr.tuned
svrPrediction.tunedRMSE <-sqrt(mean(binding.modelsvr.error.tuned^2))
print (svrPrediction.tunedRMSE)#
svrPrediction.tunedMAE <-mean(abs(binding.predict.svr.tuned-binding.train$standard_value))
print (svrPrediction.tunedMAE)#


################################BAYES############################################
samp <- sample(nrow(binding.data), 0.6 * nrow(binding.data))
binding.train <- binding.data[samp, ]#744
binding.test <- binding.data[-samp, ]#497

binding.train.net<-binding.train%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,rtb,mw_monoisotopic, aromatic_rings,full_mwt,standard_value) # 76%
binding.test.net<-binding.test%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,rtb,mw_monoisotopic, aromatic_rings,full_mwt,standard_value) # 76%


binding.hc.net <- hc(binding.train.net, debug=TRUE)#hill clinbing

# plot network
plot(binding.hc.net)

# fit the network with the data itself
binding.hc.fit <- bn.fit(binding.hc.net, binding.train.net)

binding.hc.predict <-predict(binding.hc.fit,"standard_value",binding.test.net)

plot(binding.hc.predict)
points(binding.test.net$standard_value,col="Red")

#error
binding.modelhc.error <- binding.test.net$standard_value-binding.hc.predict
binding.hc.PredictionRMSE <-sqrt(mean(binding.modelhc.error^2))
print (binding.hc.PredictionRMSE)#2.53
binding.hc.PredictionMAE <-mean(abs(binding.hc.predict-binding.test.net$standard_value))
print (binding.hc.PredictionMAE)#2.53


cbind(predicted=binding.hc.predict,actual=binding.test.net$standard_value)

###########################Multiple Linear regression##############################
plot(log(binding.train$standard_value))

binding.lin.reg <- lm(log(standard_value+1) ~ mw_freebase +  psa + acd_most_apka + acd_most_bpka + hba +
                          hbd + rtb+mw_monoisotopic+aromatic_rings+full_mwt, data = binding.train)


# Inspect the model
summary(binding.lin.reg)
#r-squared 14% ou variancia explicada pelos preditores
#PR(>|t|) é o p e por isso quanto menor mlhr (especialmente menor que 0,05)

binding.pred.lin <- exp(predict(binding.lin.reg,binding.test))-1

plot(binding.pred.lin)
points(binding.test$standard_value,col='blue')
plot(binding.pred.lin-binding.test$standard_value)
str(binding.pred.lin)

binding.lin.reg.RMSE <- sqrt(mean((binding.pred.lin-binding.test$standard_value)^2))
print(binding.lin.reg.RMSE)
binding.lin.reg.MAE <- mean(abs(binding.pred.lin-binding.test$standard_value))
print(binding.lin.reg.MAE)

################################# Neural Networks ###################################

binding.data.NN<-binding.data%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,rtb,mw_monoisotopic, aromatic_rings,full_mwt,standard_value) # 76%

#check if no na is available
apply(binding.data.NN,2,function(x) sum(is.na(x)))#fun?

binding.lm.fit<-glm(standard_value~.,data=binding.train.net)

summary(binding.lm.fit)

binding.pr.lm<-predict(binding.lm.fit,binding.test.net)
binding.MSE.lm<-sum((binding.pr.lm-binding.test.net$standard_value)^2)/nrow(binding.test)
print(binding.MSE.lm)



binding.maxs<-apply(binding.data.NN,2,max)
binding.mins<-apply(binding.data.NN,2,min)

binding.scaled<-as.data.frame(scale(binding.data.NN,center=binding.mins,scale=binding.maxs-binding.mins))

binding.train.NN<-binding.scaled[samp,]
binding.test.NN<-binding.scaled[-samp,]

binding.n<-names(binding.train.NN)
binding.f<-as.formula(paste("standard_value ~",paste(binding.n[!binding.n %in% "standard_value"],collapse = " + ")))
#o y~ n funciona no neuralnet, dai se fazer a formular fora

binding.NN.model<-neuralnet(binding.f,data=binding.train.NN,hidden=c(6,3),linear.output = TRUE)
#algorithm did not converge in 1 of 1 repetition(s) within the stepmax 
#aumentar o stepmax para dar tempo a NN de 

#You can increase the stepmax and thereby giving it more time to converge.
#The other option is to adjust the threshold parameter. 
#By default its value is 0.01. Try increasing it to 0.1/0.5. 
#If you change the lifesign to 'full' you can see the threshold values. 
#Keep your threshold value lower than the one you see at the last step. 
#Remember, higher the threshold, lower the accuracy of the model

plot(binding.NN.model)

str(binding.test.NN[,1:10])#retirar o standard_value
binding.pred.NN<-compute(binding.NN.model,binding.test.NN[,1:10])

binding.pred.NN_<-binding.pred.NN$net.result*(max(binding.data.NN$standard_value)-min(binding.data.NN$standard_value))+min(binding.data.NN$standard_value)
binding.pred.NN.response<-(binding.test.NN$standard_value)*(max(binding.data.NN$standard_value)-min(binding.data.NN$standard_value))+min(binding.data.NN$standard_value)



binding.nn.RMSE<-sum((binding.pred.NN.response-binding.pred.NN_)^2)/nrow(binding.test.NN)
print(binding.nn.RMSE)

plot(binding.test.NN$standard_value,binding.pred.NN$net.result,col="blue")
abline(0,1,lwd=2)


#function for cross-validation Neural Network.
binding.cv.NN.error<-NULL
k<-10

pbar <- create_progress_bar('text')
pbar$init(k)

for(i in 1:k){
  index <- sample(1:nrow(binding.data),round(0.9*nrow(binding.data)))
  binding.train.NN.cv <- binding.scaled[index,]
  binding.test.NN.cv <- binding.scaled[-index,]
  
  binding.nn.cv <- neuralnet(binding.f,data=binding.train.NN.cv,hidden=c(6,2),linear.output=T)
  
  binding.pr.nn.cv <- compute(binding.nn.cv,binding.test.NN.cv[,1:10])
  binding.pr.nn.cv <- binding.pr.nn.cv$net.result*(max(binding.data$standard_value)-min(binding.data$standard_value))+min(binding.data$standard_value)
  
  binding.test.cv.r <- (binding.test.NN.cv$standard_value)*(max(binding.data$standard_value)-min(binding.data$standard_value))+min(binding.data$standard_value)
  
  binding.cv.NN.error[i] <- sum((binding.test.cv.r - binding.pr.nn.cv)^2)/nrow(binding.test.NN.cv)
  
  pbar$step()
}

mean(binding.cv.NN.error)

boxplot(binding.cv.NN.error,xlab='MSE NN CV',col='green',
        border='blue',names='CV error (MSE)',
        main='CV error (MSE) for NN',horizontal=TRUE)


######################################TOTAL#######################################
binding.accuracy <- data.frame(Method = c("Random Forest","SVR","HC.Bayes","SVR Tuned","Random forest Tuned","Multiple Linear Regression","Neural Networks"),
                       RMSE   = c(rf.predction.RMSE,binding.svrPredictionRMSE,binding.hc.PredictionRMSE,svrPrediction.tunedRMSE,binding.rf.Prediction.tunedRMSE,binding.lin.reg.RMSE,mean(binding.cv.NN.error)),
                       MAE    = c(rf.predction.MAE,binding.svrPredictionMAE,binding.hc.PredictionMAE,svrPrediction.tunedMAE,binding.rf.Prediction.tunedMAE,binding.lin.reg.MAE,0)) 



# Round the values and print the table
binding.accuracy$RMSE <- round(binding.accuracy$RMSE,2)
binding.accuracy$MAE <- round(binding.accuracy$MAE,2) 

print (binding.accuracy)

######################################OUTROS#######################################

###rf Cross validation 2
binding.rf.cv <- rf.crossValidation(binding.rf.mdl, binding.x.train, p=0.60, n=10, ntree=1500)# segunda versao cv
print(binding.rf.cv)

binding.rf.cv.tuned <- rf.crossValidation(binding.rf.mdl.tuned, binding.x.train, p=0.10, n=10, ntree=1500)# segunda versao cv
print(binding.rf.cv.tuned)

