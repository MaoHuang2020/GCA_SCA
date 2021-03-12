#Run_between_loc_yBLUE_as_y
rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal

load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot

Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])


### Get the "BLUE_Trait"
load(paste0(WD,"OneTime1920/data/","BLUE_DwPM_2vs1Year_Update03082021.rdata"))

##!!!!!!!!!! BLUE is estimated within 2019 and 2020 !

for (Yr in c(2019)) {
  Y1<-Y[Y$Year==Yr,]
  Y1$BLUE_Trait<-expss::vlookup(Y1$Crosses,dict=CrossBLUE,result_column = "BLUE_DwPM_2019",lookup_column = "CrossName")  # This is the DwPM
} 

for (Yr in c(2020)){
  Y2<-Y[Y$Year==Yr,]
  Y2$BLUE_Trait<-expss::vlookup(Y2$Crosses,dict=CrossBLUE,result_column = "BLUE_DwPM_2020",lookup_column = "CrossName") 
}

Yrbind<-rbind(Y1,Y2)
Ydata<-Y
Y<-Yrbind

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


reps<-20 # !!!
ntraits<-1  # !!!
head(Y)  
y<-yBLUE<-Y[,"BLUE_Trait"] 


#CB=Farmed skinny kelp from giant staircase 
#CC: cape Cod Canal; JS:; LD:Lubec Dock; NC:Newcastle; NL: Nubble Light; OI: Orr's Island; SF:
#(NB=narrow blade=GS=Giant staircase. TI and CB are farmed)
# JS= Fort Stark, New Hampshire (SNE?)
# SF:Sullivan Falls
##From XW PCA: NL,JS,NC,IS is grouped, close to CC

locs<-c("CC","NC","NL","JS","OI","SF","LD","CB")

cor<-matrix(nrow=reps,ncol=length(locs))
for (i in 3:reps){
  setwd(paste0(WD,"OneTime1920/BetweenLoc_output/")) #where to create Rep# folder
  dir.create(paste0("yBLUEs_Rep",i))
  WDloc<-paste0(WD,"OneTime1920/BetweenLoc_output/yBLUEs_Rep",i,"/")  # the path following WD
  
  setwd(WDloc)   #where to create loc# folder
  r<-NULL
  for(loc in locs) {
    
    dir.create(paste0("loc",loc))
    savepath<-paste0("OneTime1920/BetweenLoc_output/yBLUEs_Rep",i,"/loc",loc,"/")
    
    yNA<-y 
    testing<-which(Y$femaParLoc==loc)  # 156 plots
    yNA[testing]<-NA  # Other locations to predict this testing one
    
    fm<-BGLR::BGLR(y=yNA,
                   ETA=ETA,
                   nIter=80000,
                   burnIn=60000,
                   saveAt=paste0(WD,savepath,"Yr20To19_rep",i),
                   verbose=TRUE)
    save(fm,file=paste0(WD,savepath,"fm_Yr20To19_rep",i,".rda"))
    
    yPred<-fm$ETA[[3]]$u+fm$ETA[[2]]$u+fm$ETA[[1]]$u  #SCA+GCA2+GCA1
    predict<-data.frame(testing,
                        Crosses=Y[c(testing),]$Crosses,
                        yBLUE=yBLUE[testing],
                        yPred=yPred[testing],
                        yHat=fm$yHat[testing],
                        popChk=Y[c(testing),"popChk"],
                        Year=Y[c(testing),]$Year)
    predict<-droplevels(predict)
    r<-c(r,cor(predict$yBLUE,predict$yPred,use="complete"))
    write.table(predict,file=paste0(WD,savepath,"predictions_rep",i,".csv"),row.names=FALSE,sep=",") 
  }
  
  names(r)<-locs
  write.csv(r, paste0("r_1Loc_PP_Rep",i,".csv"))
}

print(r)

rAll<-NULL
for (i in 1:reps){
  rAll[[i]]<-read.csv(paste0(WD,"OneTime1920/BetweenLoc_output/yBLUEs_Rep",i,"/","r_1Loc_PP_Rep",i,".csv"),row.names=1)
}  

rowMeans(do.call(cbind.data.frame, rAll))

