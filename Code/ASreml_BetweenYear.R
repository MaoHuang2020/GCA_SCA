###ASreml-R trait genetic correlations between years
# save(dataNHpi,All_grm,Trait_grm,file="Asreml_OtherTraits.Rdata")

###############################################  
rm(list=ls())
#WD<-"/Users/maohuang/Desktop/Kelp/GCA_SCA/" # local
#datafdr<-"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/data/"

WD<-"/local/workdir/mh865/GCA_SCA/"
datafdr<-paste0(WD,"OneTime1920/data/")

###1. Input outCovComb4_dipOrder
load(paste0(WD,"OneTime1920/data/","outCovComb_dip_0116_2021.Rdata"))
All_grm<-outCovComb4_dipOrder  # 866

###2. Input dataNH

#### Plot level
load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23_UpdateAsh_0309_2021.rdata"))  ## Plot -- Updated Ash
load(paste0(datafdr,"dataNHim_withChk_3_sets_PhotoScore0123.rdata"))  ## Individual

dataNHpi<-dataNHpiBoth_C  ##!!!!!

exptlSP <- as.character(dataNHpi$popChk) == "ES"
dataNHpi$entry <- as.character(dataNHpi$popChk)
dataNHpi$entry[exptlSP] <- as.character(dataNHpi$plotNo[exptlSP])
dataNHpi$group <- as.factor(ifelse(exptlSP, 1, 0))

for (col in c( "line", "block","popChk","group","entry","Year")){ 
  dataNHpi[,col] <- factor(dataNHpi[,col])}
dataNHpi<-droplevels(dataNHpi)

### Plot Level----Add blades!!!!
traits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades") 
dataNH<-dataNHpi   ### !!!!!


#### Individual level to get their experimental factors.
#### The dataNHpi is already filtered for their phenotypic PHOTO SCORE (>2)
dataNHim<-dataNHimboth_C
dataNHim$line<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "line")
dataNHim$block<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "block")
dataNHim$Year<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Year")
dataNHim$popChk<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "popChk")
dataNHim$Crosses<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Crosses")
dataNHim$PhotoScore<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "PhotoScore")
dataNHim<-dataNHim[which(dataNHim$PhotoScore >1),]  # 3969 rows with PhotoScore >1
  str(dataNHim)
  dim(dataNHim)

traits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")

for (col in traits){
  dataNHim[,col]<-as.numeric(dataNHim[,col]) }

dataNH<-dataNHim   ### !!!!!
#######


#R
library(asreml)
asreml.license.activate()
#enter this code CDEA-HECC-CDAH-FIED

modSum<-NULL
cor<-NULL
for (j in 1:length(traits)){

dataNH$Trait<-dataNH[,traits[j]]
dataNH<-droplevels(dataNH[!is.na(dataNH$Trait),])

if (traits[j]%in%c("dryWgtPerM","AshFDwPM")){
  dataNH$Trait2<-ifelse(dataNH$Year==2019,dataNH$Trait*sqrt(10),dataNH$Trait)  # Yr1 phenotype * sqrt(10)
  dataNH$Trait<-dataNH$Trait2
}

Trait_Crosses<-as.character(dataNH$Crosses)
Trait_grm<-All_grm[rownames(All_grm)%in%Trait_Crosses,colnames(All_grm)%in%Trait_Crosses]
    print(dim(dataNH))
    print(sum(is.na(dataNH$Trait)))
    print(dim(Trait_grm))

### Adding the checks into the Trait_grm, all 1s in diagonal,all 0s for others
data<-dataNH    
    data[data$popChk=="ES",]$Crosses
    droplevels(data[data$popChk=="ES",])$Crosses
    ChkCross<-unique(droplevels(data[!data$popChk=="ES",])$Crosses) # 33 plots of checks, 4 unique ones
    Col0<-matrix(0,nrow=nrow(Trait_grm),ncol=length(ChkCross))
    colnames(Col0)<-ChkCross
    
    Trait_grm2<-cbind(Trait_grm,Col0)
    
    Row0<-matrix(0,nrow=length(ChkCross),ncol=ncol(Trait_grm))
    Chk1<-diag(x=1,nrow=length(ChkCross),ncol=length(ChkCross))
    Row0_Chk1<-cbind(Row0,Chk1)
    rownames(Row0_Chk1)<-ChkCross
    
Trait_grm3<-rbind(Trait_grm2,Row0_Chk1)
    
  print(nrow(Trait_grm3)==length(unique(droplevels(data$Crosses))))
modSum[[j]]<-summary(Gencor_Yr(data=data,Trait_grm=Trait_grm3)$mod)$varcomp
#covarianceA_B/sqrt(varianceA*varianceB)

cor<-c(cor,modSum[[j]][,"component"][2]/sqrt(modSum[[j]][,"component"][1]*modSum[[j]][,"component"][3]))

}

names(cor)<-traits

Plotcor<-cor
write.csv(cor,"Genetic_Correlations_BetweenYears_PlotLevel.csv")
save(modSum,file="Genetic_Correlations_BetweenYears_PlotLevel_moddelSummary.Rdata")

Indicor<-cor
write.csv(cor,"Genetic_Correlations_BetweenYears_IndivLevel.csv")
save(modSum,file="Genetic_Correlations_BetweenYears_IndivLevel_modelSummary.Rdata")

Bothcor<-c(Plotcor,Indicor)
  Bothcor
write.csv(Bothcor,"Genetic_Correlations_BetweenYears_Plot_Indiv_Level.csv")

Gencor_Yr<-function (data=dataNH,Trait_grm=Trait_grm3){
mod4<-asreml(Trait~Year+line+block+popChk,
             random= ~ us(Year):vm(Crosses,Trait_grm3),
             data = data, maxiter=100, trace=TRUE)
return(list(mod=mod4))
  }

####### Multi-Trait
source("/local/workdir/mh865/GCA_SCA/OneTime1920/code/is.positive.definite.R")
source("/local/workdir/mh865/GCA_SCA/OneTime1920/code/is.square.matrix.R")
source("/local/workdir/mh865/GCA_SCA/OneTime1920/code/is.symmetric.matrix.R")
# 

# Y<-dataNHpi

#### Use the de-regressed BLUPs for multi-trait analysis?
datafdr<-paste0(WD,"OneTime1920/data/")
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear.Rdata")) ##!!!

Y<-WithinYr_Both_dBLUPs[,!colnames(WithinYr_Both_dBLUPs)%in%c("Row.names","Year.y","plotNo.y","Crosses.y")]
rownames(Y)<-Y$plotNo.x
Y$Crosses<-as.factor(Y$Crosses.x)
Y$Year<-Y$Year.x

is.positive.definite(outCovComb4_dipOrder)
outCovComb4_dipOrder[1:4,1:5]

Amat<-outCovComb4_dipOrder[rownames(outCovComb4_dipOrder)%in%as.character(Y$Crosses),colnames(outCovComb4_dipOrder)%in%as.character(Y$Crosses)]

snpRelMat<-Amat
Gsnp=solve(snpRelMat+diag(1e-6,length(snpRelMat[,1]),length(snpRelMat[,1])))

# Map elements in the relationship matrix to the phenotypes
  nrow(Gsnp)==length(unique(Y$Crosses))
  sum(rownames(Gsnp)==levels(Y$Crosses))
  sum(colnames(Gsnp)==levels(Y$Crosses))

Gsnp2<-Gsnp[match(levels(Y$Crosses),rownames(Gsnp)),match(levels(Y$Crosses),colnames(Gsnp))]
Gsnp<-Gsnp2

attr(Gsnp, "INVERSE")=TRUE

for (t in c("dryWgtPerM","AshFDwPM")){
  Y$Trait<-ifelse(Y$Year==2019,Y[,t]*sqrt(10),Y[,t])  # Yr1 phenotype * sqrt(10)
  colnames(Y)[colnames(Y)=="Trait"]<-paste0(t,"_sqrt10")
}

# running two different GxE models
modMTM <- asreml(cbind(bladeLength,bladeMaxWidth,bladeThickness,stipeLength,stipeDiameter) ~ trait, 
                 random= ~ us(trait):vm(Crosses,Gsnp),
                 residual = ~id(units):us(trait),
                 data = Y, maxiter=50, trace=TRUE)

modMTM2 <- asreml(cbind(dryWgtPerM,bladeLength,bladeMaxWidth,bladeThickness,stipeLength,stipeDiameter) ~ trait, 
                 random= ~ us(trait):vm(Crosses,Gsnp),
                 residual = ~id(units):us(trait),
                 data = Y, maxiter=50, trace=TRUE)

#generating estimates of Crosses performance in each trait
predMTM = predict(modMTM, classify = "trait:Crosses", trace=F)

# pulling out the predictions and comparing to the true simulated values
pMTM=predMTM$pvals$predicted.value

GenVar<-summary(modMTM)$varcomp
write.csv(GenVar,"Individual_MultiTrait_Genetic_var.csv")

getwd()

##Manually replaced extra char strings
Indi_Genetic_Cor<-read.csv(paste0(datafdr,"Individual_MultiTrait_Genetic_var.csv"),row.names=1)
head(Indi_Genetic_Cor)

grepl(":",Indi_Genetic_Cor$X.1)
Indi_traits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")
rownames(Indi_Genetic_Cor)[1:(3*length(Indi_traits))]
  #coVarA_B/sqrt(VarA*VarB)

for (i in 1:length(Indi_traits)){
  trait<-Indi_traits[i]
  traitname<-Indi_Genetic_Cor$X.1
  Itself[i]<-Indi_Genetic_Cor[,"component"][grepl(paste0(trait,":",trait),traitname)]
  WithOther[i]<-Indi_Genetic_Cor[,"component"][grepl(paste0(trait,":",))]
}