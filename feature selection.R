#feature selection
library(caret)
print(names(data))
# load the data
str(feature.data)

feature.all<-data%>%select(structure_type,molecule_type,ro3_pass,full_mwt,molecular_species,mw_freebase,psa,acd_most_apka,full_mwt,mw_freebase,psa,acd_most_apka,
                            hba,hbd,acd_logp,acd_logd,rtb,acd_most_bpka,mw_monoisotopic,aromatic_rings,
                            confidence_score,hba_lipinski,hbd_lipinski,
                            heavy_atoms,standard_flag,num_lipinski_ro5_violations,
                            qed_weighted,alogp,parenteral,topical,oral,chirality,prodrug,natural_product,assay_id,assay_id.1,standard_value,description)
#processar os dados
preProcValues <- preProcess(feature.data, method = c("knnImpute","center","scale"))
preProcValues <- preProcess(feature.data, method = c("center","scale"))

feature.all <- data.frame(predict(dmy, newdata = feature.all))
print(preProcValues)
library('RANN')
train_processed <- predict(preProcValues, feature.data)

dmy <- dummyVars(" ~ .", data = feature.data,fullRank = T)
print (dmy)
feature.data <- data.frame(predict(dmy, newdata = feature.data))

###########
feature.data<-feature.all%>%filter(assay_id==16590|assay_id==20401)#logp

feature.data<-feature.all%>%filter(description=='Aqueous solubility')#aqueous solubility

feature.data<-feature.all%>%filter(assay_id.1==225074)#Binding energy by using the equation deltaG obsd = -RT ln KD

feature.data<-feature.data%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,acd_logp,acd_logd,heavy_atoms,rtb,hba_lipinski,mw_monoisotopic, aromatic_rings,full_mwt,mw_monoisotopic,hbd_lipinski, standard_value) 

feature.data<-feature.data%>%select(structure_type,molecule_type,ro3_pass,full_mwt,molecular_species,mw_freebase,psa,acd_most_apka,
                           hba,hbd,acd_logp,acd_logd,rtb,acd_most_bpka,mw_monoisotopic,aromatic_rings,
                           hba_lipinski,hbd_lipinski,
                           heavy_atoms,num_lipinski_ro5_violations,
                           qed_weighted,alogp,parenteral,topical,oral,chirality,prodrug,natural_product,standard_value)

str(train_transformed)



feature.data<-feature.data%>%filter(!is.na(standard_value))
feature.data<-feature.data%>%filter(!is.na(mw_freebase))
feature.data<-feature.data%>%filter(!is.na(psa))
feature.data<-feature.data%>%filter(!is.na(acd_most_apka))
feature.data<-feature.data%>%filter(!is.na(hba))
feature.data<-feature.data%>%filter(!is.na(hbd))
feature.data<-feature.data%>%filter(!is.na(acd_logp))
feature.data<-feature.data%>%filter(!is.na(acd_logd))
feature.data<-feature.data%>%filter(!is.na(rtb))
feature.data<-feature.data%>%filter(!is.na(acd_most_bpka))
feature.data<-feature.data%>%filter(!is.na(mw_monoisotopic))
feature.data<-feature.data%>%filter(!is.na(aromatic_rings))
feature.data<-feature.data%>%filter(!is.na(full_mwt))


# calculate correlation matrix
correlationMatrix <- cor(feature.data[,1:16])
# summarize the correlation matrix
print(feature.data)
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
# print indexes of highly correlated attributes
print(highlyCorrelated)

# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(standard_value~., data=feature.data, method="svmRadial", preProcess="scale", trControl=control, na.action = na.omit )
#svmLinear2
model <- train(standard_value~., data=train_transformed, method="svmRadial", trControl=control, na.action = na.omit)


# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)

control.data<-logp.data%>%select(mw_freebase, psa,acd_most_apka, acd_most_bpka, hba,hbd,rtb,mw_monoisotopic, aromatic_rings,full_mwt,standard_value) # 76%
control.data<-feature.data


# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(control.data[,1:37], control.data[,38], sizes=c(1:37), rfeControl=control )
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))



###############################tune caret GBM#####
fitControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5)

modelLookup(model='gbm')
modelLookup(model='rf')
modelLookup(model='svmRadialCost')

grid <- expand.grid(n.trees=c(10,20,50,100,500,1000),shrinkage=c(0.01,0.05,0.1,0.5),n.minobsinnode = c(3,5,10),interaction.depth=c(1,5,10))

model_gbm<-train(control.data[,1:10], control.data[,11],method='gbm',trControl=fitControl,tuneGrid=grid)

print(model_gbm)

plot(model_gbm)

model_gbm<-train(control.data[,1:10], control.data[,11],method='gbm',trControl=fitControl,tuneLength=10)
print(model_gbm)



##########Tune caret RF ########
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="random")
mtry <- sqrt(ncol(control.data))
metric <- "RMSE"
rf_random <- train(control.data[,1:37], control.data[,38], method="rf", metric=metric, tuneLength=15, trControl=control)
print(rf_random)
plot(rf_random)
# Random Search
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="random")
set.seed(seed)
mtry <- sqrt(ncol(control.data))
rf_random <- train(control.data[,1:37], control.data[,38], method="rf", metric=metric, tuneLength=15, trControl=control, ntree=1500)
print(rf_random)
plot(rf_random)
#grid search
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(2:25))
rf_gridsearch <- train(control.data[,1:37], control.data[,38], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch)
