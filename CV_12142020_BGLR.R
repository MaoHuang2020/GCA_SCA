### Make EVD data

### parameters.R Script
setwd('/local/workdir/mh865/GCA_SCA/CVData1920/output')
phenotype.file <- '../data/CrossMerge1920.csv'          # path to pehnotype file
mm.file        <- '../data/A.csv'          # path to covariates file
weight.file     <- NULL    # path to weight file if NULL no weights are used


ctr 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
std 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
weighting 	<- !is.null(weight.file) 	# Valid only if mm.file is not NULL



colIDy <- 4  # !!!!! and change to 4

#the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
# or NULL if ID is not used
# 4 #(basicallythe P2 column ---- 1126Notes)


### Making the EVD for P1 and P2
### Checks parameters
# cat('=1===> Reading and Checking Parameters');cat('\n')
# source('../input/parameters.R') 
# if(is.null(mm.file)){ stop('Please provide path and filename of molecular marker file') }
# if(is.null(mm.file)){ cat('Are you working with a grouping factor?. If not Please provide path and filename of covariate matrix file') }


### Reads data and perform consistency checks
cat('=2===> Reading Data ');cat('\n')
if(!is.null(colIDy)){
  Y <- read.table(phenotype.file,sep=',', header=TRUE, stringsAsFactors=FALSE )
  IDs<-Y[,colIDy] # All IDs need to be in the G matrix
  if(is.null(mm.file))Y=Y else rm(Y)
}

G <- read.table(mm.file,sep=',',header=TRUE,stringsAsFactors=FALSE)
name <- G[,1]
G <- G[,-1]
colnames(G) <- name
rownames(G) <- name
G <- 2*as.matrix(G) 
n<- dim(G)[1]

Index<-rownames(G)%in%IDs  ## Adding index to ensure G has the list of individuals in the IDs (P1 list)
G<-G[Index,Index]  # subset the G to match up list of P1 names

# Set up the paramter for P1 and P2 seperately and the main code too

if(!is.null(colIDy)){ stopifnot(all(IDs%in%rownames(G))) }

if(!is.null(colIDy)){
  IDs<-factor(IDs,levels=rownames(G))
  Z<-as.matrix(model.matrix(~IDs-1))  # expand the G matrix into Gc_P1 based
  G<-tcrossprod(tcrossprod(Z,G),Z)  # ZGZ'
}


save(G,file='../GP2/G.rda')  # GP2/

EVD<-eigen(G)
rownames(EVD$vectors)<-rownames(G)
## If you provide G, the RHKS will internally compute its eigen values
## provide to save time

save(EVD,file='../GP2/EVD.rda') # GP2/


######## Generating the EVD for P1xP2

G1.file <- '../GP1/G.rda'          # path to matrix file 1
G2.file <- '../GP2/G.rda'          # path to matrix file 2


## Parameters
# source('../input/parameters.R')

## Checks parameters
# cat('=1===> Reading and Checking Parameters');cat('\n')
# source('../input/parameters.R') 
#   if(is.null(G1.file) | is.null(G2.file)){ cat('Please provide path and filename G matrices 1 and 2') }
#   cat('The Eigen value decomposition will be perfomed with components \n')


load(G1.file)
G1 <- G
load(G2.file)
G2 <- G

GI <- G1*G2
EVD <- eigen(GI)

save(GI,file='../GP1P2/G.rda')
save(EVD,file='../GP1P2/EVD.rda')



### CV scheme parameters.R
setwd('/local/workdir/mh865/GCA_SCA/CVData1920/output')

folds   <- 1:10    ### If no CV, then set a -999 for "folds" in parameters
nIter  <- 50000
burnIn <- 40000
phenotype.file <- '../data/CrossMerge1920.csv' 
AB <- list() 
AB[[1]] <- '../GP1/EVD.rda'       # path to pehnotype file 
AB[[2]] <- '../GP2/EVD.rda'       # path to pehnotype file 
AB[[3]] <- '../GP1P2/EVD.rda'       # path to pehnotype file
#AB[[4]] <- '../../../../../../G/E/output/EVD.rda'       # path to pehnotype file

type <- c('RKHS','RKHS','RKHS') #!!!
colENV  <- NULL  # column in phenotype file that gives the id of the environment 
colVAR  <- 2  # column in phenotype file that gives the id of the variety 
colPhen <- 5  # phenotypes column
colCV <- 7   # CV column
CV0 <- FALSE 
ESC <- FALSE 
r <- 1
set.seed(1)


#### CV scheme main code

#Load the BGLR library
library(BGLR)

## Parameters
#source('../input/parameters.R')

Y <- read.csv(phenotype.file,sep=',', header=TRUE, stringsAsFactors=FALSE )
rownames(Y)<-Y$GID
	
#load(files)

y   = Y[,colPhen]
gid = Y[,colVAR]

if(ESC) { y=scale(y,center=TRUE,scale=TRUE) }


## This is to define what model do you want to use, in BGLR, setting up the ETA respectively
###########
#### WHY??? This is to allow each AB a separate ETA?

nk <-length(AB)  # what is the nk here? AB has GP1, GP2, GP1P2
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
########

### Setting the CV
########



#### make more random samples
sampleCV<-matrix(nrow=nrow(Y),ncol=500)
for (n in 1:500){
sets<-rep(1:10,26)[-c(1:4)]  # nrow(Y)=nrow(CrossMerge1920)=256 lines !!!!!!!!
sampleCV[,n]<-sets[order(runif(nrow(Y)))]
}

save(sampleCV,file="sampleCV.Rdata")

cycles<-50 # !!!
ntraits<-1  # !!!
cor<-matrix(nrow=cycles,ncol=ntraits)

for (i in 1:cycles){

	tmp<-NULL

for(fold in folds){

    yNA<-y
   # print(fold)


    colCV<-i
    testing=which(sampleCV[,colCV]==fold)
    
   
 #   if(CV0)
 #   {         
 #     testing <- which(gid %in% gid[testing])   #CV0 is DJ paper, random remove one at a time
 #   }  
   
  

    yNA[testing]=NA
    
    fm=BGLR(y=yNA,
	ETA=ETA,
	nIter=nIter,
	burnIn=burnIn,
	saveAt=paste0("CVData1920_10fold_Cycle",i),
	verbose=TRUE)

    	fm$y=y     ### !!!! Not clear about this step
    
    predictions=data.frame(testing,Individual=gid[testing], y=y[testing], yHat=fm$yHat[testing])
	tmp<-rbind(tmp,predictions)
}  # End folds
	cor[i,]<-cor(tmp$y,tmp$yHat)

### Save the predictions from full 10 folds
   
dir.create(paste0('10folds_Cycle',i))
setwd(paste0('10folds_Cycle',i))   ### Creating the fold_# folder
write.table(tmp,file=paste("predictions_cycle",i,".csv",sep=""),row.names=FALSE,sep=",") 
             
  # print(str(fm))
  rm(fm) 
  
  unlink("*.dat")
  
  setwd('..')

} # End cycles