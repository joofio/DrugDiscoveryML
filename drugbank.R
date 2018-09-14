#drugbank teste WIP


dbdata<-read.csv(file="drugbank.csv",sep=";")

str(dbdata)
glimpse(dbdata)

dbdata$published_value <-as.numeric(as.character(dbdata$published_value))
dbdata$Water.Solubility.EXP <-as.character(dbdata$Water.Solubility.EXP)

dbdata<-dbdata%>%select('drugbank_id', 'name', 'type',  'average.mass', 'monoisotopic.mass', 'volume.of.distribution', 'clearance', 'half.life', 'toxicity', 'metabolism', 'metabolism.1', 'absorption', 'smiles', 'melting.point', 'logS.EXP', 'logP.EXP', 'pKa.EXP', 'Isoelectric.Point', 'Molecular.Weight', 'Molecular.Formula', 'Hydrophobicity', 'Boiling.Point', 'Caco2.Permeability', 'Water.Solubility.EXP', 'PSA.calc', 'Refractivity.calc', 'Polarizability', 'Ghose.Filter', 'MDDR.Like.Rule')

dbcolumns<-names(dbdata)
print(paste(dbcolumns, collapse="', '"))

dbdata$Boiling.Point <-as.numeric(as.character(dbdata$Boiling.Point))
dbdata$name <-as.character(dbdata$name)
dbdata$Molecular.Formula <-as.character(dbdata$Molecular.Formula)
