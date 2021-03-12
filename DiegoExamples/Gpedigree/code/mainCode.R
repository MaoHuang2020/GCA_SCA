#############################################
 rm(list=ls())

## The working directory should be the analysis directory
 setwd('../output')
 setwd('C:/CIMMYT/SUGAR_CANE/ANALYSIS/TCH/Gmatrix/Gpedigree/output')
## Parameters
 source('../input/parameters.R')

## Checks parameters
 cat('=1===> Reading and Checking Parameters');cat('\n')
 source('../input/parameters.R') 
 #if(is.null(mm.file)){ stop('Please provide path and filename of molecular marker file') }
  if(is.null(mm.file)){ cat('Are you working with a grouping factor?. If not Please provide path and filename of covariate matrix file') }

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

   
   save(G,file='G.rda')

   EVD<-eigen(G)
   rownames(EVD$vectors)<-rownames(G)
## If you provide G, the RHKS will internally compute its eigen values
## provide to save time

   save(EVD,file='EVD.rda')










 pdf(file='Eigen.pdf')
  plot(EVD$vectors[,1],EVD$vectors[,2],main='',xlab='First Component',ylab='Second Component')
  plot(EVD$vectors[,1],EVD$vectors[,3],main='',xlab='First Component',ylab='Third Component')
  plot(EVD$vectors[,2],EVD$vectors[,3],main='',xlab='Second Component',ylab='Third Component')
  plot(EVD$values, main='',xlab='Components',ylab='Eigen-value')
  plot(cumsum(EVD$values)[1:I(min(n,100))]/sum(EVD$values),main='',xlab='components',ylab='Cumulative Variance')
  abline(h=.8,lty=3)
  abline(v=which(  floor(cumsum(EVD$values)[1:100]/sum(EVD$values)*10) ==8)[1], lty=3)
  dev.off()

quit(save='no')


