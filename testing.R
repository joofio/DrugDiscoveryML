#processar os dados
#preProcValues <- preProcess(feature.data, method = c("knnImpute","center","scale"))
#preProcValues <- preProcess(feature.data, method = c("center","scale"))


#str(head(data))
#data%>%group_by(alogp)%>%summarise(count=n())%>%arrange(desc(count))#12323
#data%>%group_by(standard_value)%>%summarise(count=n())%>%arrange(desc(count))#19850

######deal with Nas

#data%>%distinct(alogp)%>%arrange(desc(alogp))
#colSums(is.na(data))

#variable to measure
aqsolubil.x.train<-aqsolubil.train%>%select(full_mwt,mw_freebase,psa,acd_most_apka,hba,hbd,acd_logp,acd_logd,rtb,acd_most_bpka,mw_monoisotopic,aromatic_rings) # 68%


aqsolubil.data%>%group_by(standard_units)%>%summarise(count=n())%>%arrange(desc(count))
aqsolubil.data%>%group_by(published_units)%>%summarise(count=n())%>%arrange(desc(published_units))


print(preProcValues)
train_processed <- predict(preProcValues, feature.data)
str(train_transformed)

#data%>%group_by(assay_id,description)%>%summarise(count=n())%>%arrange(desc(count))

dmy <- dummyVars(" ~ .", data = feature.data,fullRank = T)
print (dmy)
feature.data <- data.frame(predict(dmy, newdata = feature.data))

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

#svmLinear2
model <- train(standard_value~., data=train_transformed, method="svmRadial", trControl=control, na.action = na.omit)
