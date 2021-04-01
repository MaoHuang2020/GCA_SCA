### Use all data to run onetime GCA+SCA model, Get GCA and SCA var_component

##### GCA, SCA comp

rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal

load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot

Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])

Inputfiledir<-c("OneTime1920/GP1/250Individual/","OneTime1920/GP2/250Individual/","OneTime1920/GP1P2/250Individual/") 
  head(Y)
  dim(Y)

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

load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_overTwoYears.Rdata")) ##!!!
rownames(Both_dBLUPs)<-Both_dBLUPs$Row.names
Both_dBLUPs<-Both_dBLUPs[,-1]
traits<-colnames(Both_dBLUPs)
  print(traits)
CrossBLUE<-Both_dBLUPs ##!!!
##### !!!

reps<-20
cor<-matrix(nrow=reps,ncol=length(traits))
colnames(cor)<-traits

for (j in 1:length(traits)){
  Coltrait<-traits[j]
  Y$BLUE_Trait<-expss::vlookup(Y$Crosses,dict=CrossBLUE,result_column = paste0(Coltrait),lookup_column = "row.names") ### This is the DwPM
  head(Y)
  
  head(Y)
  dim(Y)
  
  Y<-Y   # Between Year !!!!
  y<-yBLUE<-Y[,"BLUE_Trait"]  # phenotypes column  !!!

#tmp<-NULL 

for (i in 1:reps){
  setwd(paste0(WD,"OneTime1920/Alldata_output/"))
  dir.create(paste0(Coltrait,"_OneTimeAll_Rep",i))
  savepath<-paste0("OneTime1920/Alldata_output/",Coltrait,"_OneTimeAll_Rep",i,"/")  # the path following WD
  
  y<-y  
  testing<-which(!is.na(Y$crossID))  # all lines
  str(testing)
  yNA<-y[testing]
  fm<-BGLR(y=yNA,
           ETA=ETA,
           nIter=80000,
           burnIn=60000,
           saveAt=paste0(WD,savepath,"OnetimeAll_rep",i),
           verbose=TRUE)
  save(fm,file=paste0(WD,savepath,"fm.rda"))
  
  yPred<-fm$ETA[[3]]$u+fm$ETA[[2]]$u+fm$ETA[[1]]$u  #SCA+GCA2+GCA1
  predict<-data.frame(testing,
                      Crosses=Y[c(testing),]$Crosses,
                      yBLUE=yBLUE[testing],
                      yPred=yPred[testing],
                      yHat=fm$yHat[testing],
                      popChk=Y[c(testing),"popChk"],
                      Year=Y[c(testing),]$Year)
  predict<-droplevels(predict)
  
  #tmp<-rbind(tmp,predict)
   
cor[i,j]<-cor(predict$yBLUE,predict$yPred,use="complete") 
  }
write.table(predict,file=paste0(WD,savepath,"Onetimepredictions_rep",i,".csv"),row.names=FALSE,sep=",")
}

colMean(cor) 


reps<-20
GCASCA<-matrix(nrow=reps,ncol=4)
for (j in 1:length(traits)){
  Coltrait<-traits[j]
varfm<-NULL #showed 5 varcomp
for (i in 1:reps){
  varcomp<-yHatVarMean_noloc(filedir = paste0("OneTime1920/Alldata_output/",Coltrait,"_OneTimeAll_Rep",i,"/"),vGCA1=1,vGCA2=2,vSCA=3,filename=paste0("OnetimeAll_rep",i))$varMean[2,]
  varfm<-rbind(varfm,varcomp)

  #varcomp2<-yHatVarMean_withloc(filedir = "OneTime1920/Alldata_output/withFMLocRep1/",vGCA1=1,vGCA2=2,vSCA=3,filename="Model_withFMLoc")
}
  GCASCA[j,]<-colMeans(varfm)
}
colnames(GCASCA)<-colnames(varfm)
rownames(GCASCA)<-traits[1:j]
save(GCASCA,paste0(WD,"OneTime1920/Alldata_output/","OnetimeAll_varcomp.csv"))



#### 4.FUNCTION 
yHatVarMean_noloc<-function(filedir,vGCA1=1,vGCA2=2,vSCA=3,filename="OnetimeAll_rep1"){
  load(paste0(WD,filedir,"fm.rda")) 
  fm<-fm
  varfm<-c(fm$varE,fm$ETA[[vGCA1]]$varU,fm$ETA[[vGCA2]]$varU,fm$ETA[[vSCA]]$varU) # varcomp output by fm
  
  varE<-scan(paste0(WD,filedir,filename,"varE.dat"))
  # varB<-scan(paste0(WD,filedir,filename,"ETA_",vB,"_varB.dat"))
  varU1<-scan(paste0(WD,filedir,filename,"ETA_",vGCA1,"_varU.dat"))
  varU2<-scan(paste0(WD,filedir,filename,"ETA_",vGCA2,"_varU.dat"))
  varU3<-scan(paste0(WD,filedir,filename,"ETA_",vSCA,"_varU.dat"))
  
  # plot(varU3,type='o',col=2,cex=.5)
  #### Save the plot!!!!!
  pdf(paste0(WD,filedir,"varE_varGCA1_2_SCA_OnetimeAll.pdf"))
  par(mfrow=c(2,2))
  plot(varE,type='o',col=2,cex=.5)
  plot(varU1,type='o',col=1,cex=.5)
  plot(varU2,type='o',col=1,cex=.5)
  plot(varU3,type='o',col=1,cex=.5)
  dev.off()
  
  varMean<-c(mean(varE),mean(varU1),mean(varU2),mean(varU3)) #calculated varcomp by hand
  names(varMean)<-c("varE","varGCA1","varGCA2","varSCA")
  
  varMean<-rbind(varMean,varfm)     
  write.csv(varMean,paste0(WD,filedir,filename,"_varcomp.csv"))
  return(list(varMean=varMean))
}