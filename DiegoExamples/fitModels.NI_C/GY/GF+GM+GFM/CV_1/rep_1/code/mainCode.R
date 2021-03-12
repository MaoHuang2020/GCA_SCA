  #Clean workspace

  rm(list=ls())

  #Set working directory
  setwd("../output/")


  #Load the BGLR library
  library(BGLR)

  ## Parameters
  source('../input/parameters.R')

  Y <- read.csv(phenotype.file,sep=',', header=TRUE, stringsAsFactors=FALSE )


  #load(files)

  y   = Y[,colPhen]
  gid = Y[,colVAR]


  if(ESC) { y=scale(y,center=TRUE,scale=TRUE) }





  nk <-length(AB)
  ETA<-list(nk)
   for(i in 1:nk){

        if(type[i]=='BRR')
          {
           load(AB[[i]])
           ETA[[i]] <- list(X=Z,model='BRR')
           rm(Z)
          }

        if(type[i]=='FIXED')
          {
           load(AB[[i]])
           ETA[[i]] <- list(X=Z,model='FIXED')
           rm(Z)
          }

        if(type[i]=='RKHS')
          {
           load(AB[[i]])
           ETA[[i]] <- list(V=EVD$vectors,d=EVD$values,model='RKHS')
           rm(EVD)
          }

        if(type[i]=='BayesA')
          {
           load(AB[[i]])
           ETA[[i]] <- list(X=X,model='BayesA')
           rm(X)
          }

        if(type[i]=='BayesB')
          {
           load(AB[[i]])
           ETA[[i]] <- list(X=X,model='BayesB')
           rm(X)
          }

        if(type[i]=='BayesC')
          {
           load(AB[[i]])
           ETA[[i]] <- list(X=X,model='BayesC')
           rm(X)
          }

        if(type[i]=='BL')
          {
           load(AB[[i]])
           ETA[[i]] <- list(X=X,model='BL')
           rm(X)
          }  
        print(i)
   }


  for(fold in folds)
     {
 
      yNA<-y
      print(fold)

      if(fold != -999)  # There is no CVs, predict all lines
         {

          dir.create(paste('fold_',fold,sep=''))
          setwd(paste('fold_',fold,sep=''))

          testing=which(Y[,colCV]==fold)
        
          if(CV0)
            {         
             testing <- which(gid %in% gid[testing])   
            }  
       
                  
         yNA=y
         yNA[testing]=NA


         fm=BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose=TRUE)
         fm$y=y

         predictions=data.frame(testing,Individual=gid[testing], y=y[testing], yHat=fm$yHat[testing])

         write.table(predictions,file=paste("predictions_",fold,".csv",sep=""),row.names=FALSE,sep=",") # Change to a unique name?

        }else{

            dir.create('fullData')   # when folds=-999
         setwd('fullData')

         fm=BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose=TRUE)
         save(fm,file='fm_full.RData')

        }          
  
         print(str(fm))
         rm(fm)

         
         unlink("*.dat")

         setwd('..')

     }
  





