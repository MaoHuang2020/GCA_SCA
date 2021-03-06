# Run_CV_no_Loc_using_yBLUE_as_y_Within2019

rm(list=ls())

WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal

load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
#load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/dataNHim_withChk_3_sets_PhotoScore0123.rdata")  ## Indi
####
Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])

library(BGLR)



########### Run only Once !!!!!!!!!1

load(paste0(WD,"OneTime1920/data/","outCovComb_dip_0116_2021.Rdata"))
#write.csv(outCovComb4_dipOrder,here("OneTime1920/data","A.csv"))
 mm.file<- paste0(WD,"OneTime1920/data/","A.csv")          # path to covariates file

## This is to only get the calc GCA and SCA functions
source(paste0(WD,"OneTime1920/code/","BGLR_functions_noFMLoc.R")) # !!!terminal

#############
 
 Y<-droplevels(Y[Y$Year==2019,])
 Yr<-2019
 folder<-"OneTime1920/Yr19Only/"
 nreps<-13
 subtract<-c(1:7)
 
 Y<-droplevels(Y[Y$Year==2020,])
 Yr<-2020
 folder<-"OneTime1920/Yr20Only/"
 nreps<-13
 subtract<-c(1:2)
#the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = which(colnames(Y)=="femaPar"),colNam = "P1",savefiledir = paste0(folder,"GP1/"))
CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = which(colnames(Y)=="malePar"),colNam = "P2",savefiledir = paste0(folder,"GP2/"))

CalSCA(G1.file=paste0(WD,folder,"GP1/","G.rda"),
       G2.file=paste0(WD,folder,"GP2/","G.rda"),
       savefileDir=paste0(folder,"GP1P2/"))
############# Run it only once

#### Run it once
sampleCV<-matrix(nrow=nrow(Y),ncol=500)
 for (n in 1:500){
   sets<-rep(1:10,nreps)[-subtract]  #122 ES ones
   sampleCV[,n]<-sets[order(runif(nrow(Y)))]
 }
save(sampleCV,file=paste0(WD,"OneTime1920/data/sampleCV_Yr",Yr,"_122Indiv_0312_2021.Rdata"))
 ###### One Time

############## Run only Once !!!!!!!!!



Yr<-2019  ##!!!
folder<-"OneTime1920/Yr19Only/"


Yr<-2020
folder<-"OneTime1920/Yr20Only/"

### Load Sample file
load(paste0(WD,"OneTime1920/data/sampleCV_Yr",Yr,"_122Indiv_0312_2021.Rdata"))

Inputfiledir<-c(paste0(folder,"GP1/"),paste0(folder,"GP2/"),paste0(folder,"GP1P2/")) 

load(paste0(WD,Inputfiledir[1],"EVD.rda"))  
EVD1<-EVD
  rm(EVD)
load(paste0(WD,Inputfiledir[2],"EVD.rda"))       
EVD2<-EVD
  rm(EVD)
load(paste0(WD,Inputfiledir[3],"EVD.rda"))       
EVD3<-EVD
  rm(EVD)
ETA<-list(
  list(V=EVD1$vectors,d=EVD1$values,model="RKHS"),
  list(V=EVD2$vectors,d=EVD2$values,model="RKHS"),
  list(V=EVD3$vectors,d=EVD3$values,model="RKHS")
)


#####!!!
datafdr<-paste0(WD,"OneTime1920/data/")
  ##!!! WithinYear Data
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata")) 
rownames(WithinYr_Both_dBLUPs)<-WithinYr_Both_dBLUPs$Row.names
WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,-1]

CrossBLUE<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==Yr,] ### Subset the Yr
CrossBLUE<-CrossBLUE[,!colnames(CrossBLUE)%in%c("Year.x","Year.y")] ### RM extra cols

traits<-colnames(CrossBLUE)
  print(traits)
##### !!!

sampleCV<-sampleCV
folds   <- 1:10 
reps<-20 # !!!
ntraits<- length(traits) # !!!
cor<-matrix(nrow=reps,ncol=ntraits)
colnames(cor)<-traits

for (j in 1:length(traits)){
  Coltrait<-traits[j]

  Y$BLUE_Trait<-expss::vlookup(Y$Crosses,dict=CrossBLUE,result_column = paste0(Coltrait),lookup_column = "row.names") ### !!!
  head(Y)
  
  head(Y)
  dim(Y)
  
  Y<-droplevels(Y[Y$Year==Yr,])   # Within Year !!

  y<-Y[,"BLUE_Trait"]  # phenotypes column  !!!

  yBLUE<-Y[,"BLUE_Trait"] # This is the BLUE for DwPM
  
  setwd(paste0(WD,"OneTime1920/Alldata_CV_output/"))

  for (i in 1:reps){
    setwd(paste0(WD,"OneTime1920/Alldata_CV_output/"))
    dir.create(paste0(Coltrait,"_OnlyYr",Yr,"_ydrBLUPsnoLocRep",i))
    savepath<-paste0("OneTime1920/Alldata_CV_output/",Coltrait,"_OnlyYr",Yr,"_ydrBLUPsnoLocRep",i,"/")  # the path within WD!
    
    tmp<-NULL
    for (fold in folds){
      
      #dir.create(paste0('10folds_Cycle',i,"/"))
      #savepath<-paste0('10folds_Cycle',i,"/")   ### Creating the fold_# folder
      
      yNA<-y
      # print(fold)
      testing<-which(sampleCV[,i]==fold)
      yNA[testing]<-NA
      
      fm<-BGLR(y=yNA,
               ETA=ETA,
               nIter=80000,
               burnIn=60000,
               saveAt=paste0(WD,savepath,"CVData1920_",fold,"thfold_rep",i),
               verbose=TRUE)
      save(fm,file=paste0(WD,savepath,"fm_",fold,"thfold_rep",i,".rda"))
      
      yPred<-fm$ETA[[3]]$u+fm$ETA[[2]]$u+fm$ETA[[1]]$u  #SCA+GCA2+GCA1
      predict<-data.frame(testing,
                          Crosses=Y[c(testing),]$Crosses,
                          yBLUE=yBLUE[testing],
                          yPred=yPred[testing],
                          yHat=fm$yHat[testing],
                          popChk=Y[c(testing),"popChk"],
                          Year=Y[c(testing),]$Year)
      predict<-droplevels(predict)
      
      tmp<-rbind(tmp,predict)
      # r<-c(r,predict(testing=testing,gid=gid,yBLUE=yBLUE,Y=Y,fmfiledir=savepath))
    }
    cor[i,j]<-cor(tmp$yBLUE,tmp$yPred,use="complete")
    
    write.table(tmp,file=paste0(WD,savepath,"predictions_rep",i,".csv"),row.names=FALSE,sep=",") 
    rm(fm) 
    unlink("*.dat")
  } 
  
}

### AddBD (change "length(traits)" to "AddBD)
write.csv(colMeans(cor),paste0(paste0("cor_CV_no_Loc_OnlyYr",Yr,"_ydrBLUPs_data_",length(traits),"Traits_Mean.csv")))
write.csv(cor,paste0("cor_CV_no_Loc_OnlyYr",Yr,"_ydrBLUPs_data_",length(traits),"Traits.csv"))
