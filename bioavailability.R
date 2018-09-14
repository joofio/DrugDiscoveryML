#Bioavailability
#very phew data 

set.seed(123)

#Aqueous solubility
bioavail.data<-data%>%filter(description=='Bioavailability')
bioavail.data<-data%>%filter(grepl('(Bioavailability)',description))#pouco

bioavail.data%>%group_by(description)%>%summarise(count=n())%>%arrange(desc(count))

bioavail.data<-bioavail.data%>%filter(standard_units!="NULL")

#remover NA
bioavail.data<-bioavail.data%>%filter(!is.na(mw_freebase))
bioavail.data<-bioavail.data%>%filter(!is.na(psa))
bioavail.data<-bioavail.data%>%filter(!is.na(acd_most_apka))#60
bioavail.data<-bioavail.data%>%filter(!is.na(hba))
bioavail.data<-bioavail.data%>%filter(!is.na(hbd))
bioavail.data<-bioavail.data%>%filter(!is.na(acd_logp))
bioavail.data<-bioavail.data%>%filter(!is.na(acd_logd))
bioavail.data<-bioavail.data%>%filter(!is.na(rtb))
bioavail.data<-bioavail.data%>%filter(!is.na(acd_most_bpka))#outros 60
bioavail.data<-bioavail.data%>%filter(!is.na(mw_monoisotopic))
bioavail.data<-bioavail.data%>%filter(!is.na(aromatic_rings))
bioavail.data<-bioavail.data%>%filter(!is.na(full_mwt))

bioavail.data<-bioavail.data%>%filter(!is.na(standard_value))


samp <- sample(nrow(bioavail.data), 0.7 * nrow(bioavail.data))
bioavail.train <- bioavail.data[samp, ]#744
bioavail.test <- bioavail.data[-samp, ]#497

bioavail.x.train<-bioavail.train%>%select(full_mwt,mw_freebase,psa,acd_most_apka,hba,hbd,acd_logp,acd_logd,rtb,acd_most_bpka,mw_monoisotopic,aromatic_rings) # 64%
bioavail.y.train<-bioavail.train$published_value

str(bioavail.y.train)

bioavail.rf.mdl <-randomForest(bioavail.y.train,x=bioavail.x.train,type=1)
print(bioavail.rf.mdl)
varImpPlot(bioavail.rf.mdl)
importance(bioavail.rf.mdl, scale = TRUE)

bioavail.pred<-predict(bioavail.rf.mdl,newdata=bioavail.test)

plot(bioavail.pred,bioavail.test$standard_value)
abline(0,1,col="blue")
###MSE propriamente dito, e n previsto como no print rf
mean((bioavail.pred-bioavail.test$standard_value)^2) #500+

