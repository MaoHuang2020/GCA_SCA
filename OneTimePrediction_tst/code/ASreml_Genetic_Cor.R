# ASReml

#R
library(asreml)
asreml.license.activate()
#enter this code CDEA-HECC-CDAH-FIED

setwd("/local/workdir/mh865/ASreml")
source("is.symmetric.matrix.R")
source("is.square.matrix.R")
source("is.positive.definite.R")

FileDir<-"/local/workdir/mh865/GCA_SCA/" 
load(paste0(FileDir,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
load(paste0(FileDir,"OneTime1920/data/","outCovComb_dip_0116_2021.Rdata"))

data<-droplevels(dataNHpiBoth_C)
data$Trait<-data$dryWgtPerM

data$Year<-as.factor(data$Year)
data$popChk<-as.factor(data$popChk)
data$line<-as.factor(data$line)
data$block<-as.factor(data$block)

Trait_grm<-outCovComb4_dipOrder[rownames(outCovComb4_dipOrder)%in%as.character(data$Crosses),colnames(outCovComb4_dipOrder)%in%as.character(data$Crosses)]

### Adding the checks into the grm, all 1s in diagonal,all 0s for others
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

data$Trait2<-ifelse(data$Year==2019,data$Trait*sqrt(10),data$Trait)  # Yr1 phenotype * sqrt(10)

data1<-droplevels(data[data$Year==2019,])
data2<-droplevels(data[data$Year==2020,])

data1_grm<-Trait_grm3[rownames(Trait_grm3)%in%as.character(data1$Crosses),colnames(Trait_grm3)%in%as.character(data1$Crosses)]

data2_grm<-Trait_grm3[rownames(Trait_grm3)%in%as.character(data2$Crosses),colnames(Trait_grm3)%in%as.character(data2$Crosses)]

modBoth <- asreml(Trait2 ~ Year+line+block+popChk,
               random= ~ us(Year):vm(Crosses,Trait_grm3),
               data = data, maxiter=100, trace=TRUE)

mod1 <- asreml(Trait ~ line+block+popChk,
               random= ~ vm(Crosses,data1_grm),
               data = data1, maxiter=100, trace=TRUE)

mod2 <- asreml(Trait ~ line+block+popChk,
               random= ~ vm(Crosses,data2_grm),
               data = data2, maxiter=100, trace=TRUE)

save(data,Trait_grm,Trait_grm3,file="dataNHpiBoth_C_for_ASReml.rdata")

summary(mod1)$varcomp
#covarianceA_B/sqrt(varianceA*varianceB)

### Multi-Trait BLUP
#1. Make Ainverse
ls()
source("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Code_10032020/is.positive.definite.R")
source("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Code_10032020/is.square.matrix.R")
source("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Code_10032020/is.symmetric.matrix.R")
 # mat2sparse <- function (X, rowNames = dimnames(X)[[1]])
 # {
 #   which <- (X != 0 & lower.tri(X, diag = TRUE))
 #   df <- data.frame(row = t(row(X))[t(which)], col = t(col(X))[t(which)],
 #                    val = t(X)[t(which)])
 #   if (is.null(rowNames))
 #     rowNames <- as.character(1:nrow(X))
 #   attr(df, "rowNames") <- rowNames
 #   df
 # }
 
# Trait<-"DwPM"
 
###### 
 for (t in c("dryWgtPerM","AshFDwPM")){
    dataNH$Trait<-ifelse(dataNH$Year==2019,dataNH[,t]*sqrt(10),dataNH[,t])  # Yr1 phenotype * sqrt(10)
    colnames(dataNH)[colnames(dataNH)=="Trait"]<-paste0(t,"_sqrt10")
   }
 
 
 Y<-dataNH

 is.positive.definite(outCovComb4_dipOrder)
   outCovComb4_dipOrder[1:4,1:5]

 Amat<-outCovComb4_dipOrder[rownames(outCovComb4_dipOrder)%in%as.character(Y$Crosses),colnames(outCovComb4_dipOrder)%in%as.character(Y$Crosses)]

 
 snpRelMat<-Amat
 Gsnp=solve(snpRelMat+diag(1e-6,length(snpRelMat[,1]),length(snpRelMat[,1])))
 
 # Map elements in the relationship matrix to the phenotypes
 
 #rownames(Gsnp)=levels(dataf$variety)
 #colnames(Gsnp)=levels(dataf$variety)
 attr(Gsnp, "INVERSE")=TRUE
 
 # running two different GxE models
 modMTM <- asreml(cbind(wetWgtPerM,percDryWgt,dryWgtPerM,densityBlades) ~ trait, 
                  random= ~ us(trait):vm(Crosses,Gsnp),
                  residual = ~id(units):us(trait),
                  data = Y, maxiter=50, trace=TRUE)
 


 
####### Older way of package syntax???

 #N<-nrow(Amat)
 # AHAT.inv<-solve(Amat)
 #   AHAT.inv[1:5,1:5]
 #   det(AHAT.inv)
 #AHAT.inv.sparse<-mat2sparse(AHAT.inv)  #lower diag sparse matrix
 # colnames(AHAT.inv.sparse)<-c('Row','Column','Ainverse')
 #   head(AHAT.inv.sparse)
 # write.table(AHAT.inv.sparse,file=paste0(Trait,"_Amat-Inv-sparse.txt"))
 #  #2. Make dummy pedi_file
 #  peddummy<-matrix(ncol=3,nrow=N)
 #  colnames(peddummy)<-c("Individual","Female","Male")
 #  peddummy[,1]<-rownames(Amat)
 #  peddummy[,2]<-rep(0,length=N)
 #  peddummy[,3]<-rep(0,length=N)
 #  peddummy<-as.data.frame(peddummy)
 #  rownames(peddummy)<-rownames(Amat)
 #  write.table(peddummy,file=paste(Trait,"_peddummy.txt",sep=""))
 # 
 #  gmatrix<-data.frame(AHAT.inv.sparse)
 #  ainvped<-ainverse(peddummy)
 #  attr(gmatrix,"rowNames")<-attr(ainvped,"rowNames")
 #  Y$Genot<-as.factor(Y$Crosses)  ###
# MultiTrait<-asreml(fixed=cbind(wetWgtPlot,wetWgtPerM,percDryWgt,dryWgtPerM,AshFreedryWgtPerM,densityBlades)~trait+line+block+popChk+Year,
#                    residual=~id(units):us(trait),
#                    random=~vm(Crosses,Trait_grm3),
#                    workspace=128e06,na.action=na.method(y="include"),
#                    data=data)
