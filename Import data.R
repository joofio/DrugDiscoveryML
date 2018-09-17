library(party)
library(randomForest)
library(plyr)
library(dplyr)
library(rfUtilities)
library(ggplot2)
library (e1071)
library(bnlearn)
library(neuralnet)

source("http://bioconductor.org/biocLite.R")
biocLite(c("graph","RBGL","Rgraphviz"))
install.packages(c("bnlearn","gRain","pROC","epitools"))

#set seed
set.seed(123)

#read csv origin Chembl23
data<-read.csv(file="/product_adme.csv")

str(data)

ToSelect<-data%>%group_by(assay_id,description)%>%summarise(count=n())%>%arrange(desc(count))

#correction of data types
data$standard_value <-as.numeric(as.character(data$standard_value))
data$alogp <-as.numeric(as.character(data$alogp))
data$hba <-as.numeric(as.character(data$hba))
data$hbd <-as.numeric(as.character(data$hbd))
data$mw_freebase <-as.numeric(as.character(data$mw_freebase))
data$acd_most_apka <-as.numeric(as.character(data$acd_most_apka))
data$acd_most_bpka <-as.numeric(as.character(data$acd_most_bpka))
data$num_lipinski_ro5_violations <-as.numeric(as.character(data$num_lipinski_ro5_violations))
data$psa <-as.numeric(as.character(data$psa))
data$acd_logp <-as.numeric(as.character(data$acd_logp))
data$acd_logd <-as.numeric(as.character(data$acd_logd))
data$heavy_atoms <-as.numeric(as.character(data$heavy_atoms))
data$qed_weighted <-as.numeric(as.character(data$qed_weighted))
data$num_ro5_violations <-as.numeric(as.character(data$num_ro5_violations))
data$rtb <-as.numeric(as.character(data$rtb))
data$hba_lipinski <-as.numeric(as.character(data$hba_lipinski))
data$mw_monoisotopic <-as.numeric(as.character(data$mw_monoisotopic))
data$full_mwt <-as.numeric(as.character(data$full_mwt))
data$aromatic_rings <-as.numeric(as.character(data$aromatic_rings))
data$mw_monoisotopic <-as.numeric(as.character(data$mw_monoisotopic))
data$hbd_lipinski <-as.numeric(as.character(data$hbd_lipinski))
data$published_value <-as.numeric(as.character(data$published_value))
data$standard_units <-as.character(data$standard_units)
data$description <-as.character(data$description)



test<-data%>%group_by(alogp)%>%summarise(count=n())%>%arrange(desc(count))#12323
test<-data%>%group_by(standard_value)%>%summarise(count=n())%>%arrange(desc(count))#19850


#str(head(data))

######deal with Nas

#data%>%distinct(alogp)%>%arrange(desc(alogp))
#colSums(is.na(data))

