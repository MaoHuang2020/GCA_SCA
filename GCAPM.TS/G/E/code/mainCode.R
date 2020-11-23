#############################################
 rm(list=ls())

## The working directory should be the analysis directory
 setwd('../output')

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
   IDs<-Y[,colIDy]
   if(is.null(mm.file))Y=Y else rm(Y)
 }

if(is.null(mm.file))
{
 Y[,colIDy]<-factor(Y[,colIDy])
 n<-length(levels(Y[,colIDy]))
 Z<-as.matrix(model.matrix(~Y[,colIDy]-1))
 d<-colSums(Z)
 V<-Z
 for(i in 1:ncol(Z)){ V[,i]<-V[,i]/sqrt(d[i]) } 
 EVD<-list(vectors=V,values=d)
 G <- tcrossprod(Z) 
   save(G,file='G.rda')
   save(EVD,file='EVD.rda')
}else{

   ## Reads covariates
   n<-length(count.fields(mm.file))-1


   fileIn<-file(mm.file,open='r')
   colNames<-scan(fileIn,nlines=1,what=character(),quiet=TRUE, sep=',')
   if(!is.null(colIDy)){colNames=colNames[-1]}
 
   p<-length(colNames)
 
   X<-matrix(nrow=n,ncol=p,NA)
   IDx<-rep(NA,n)
   for(i in 1:n){
      tmp<-scan(fileIn,nlines=1,what=character(),quiet=TRUE,sep=',')
      if(!is.null(colIDy))  IDx[i]<-tmp[1]
      ifelse(!is.null(colIDy)==FALSE, X[i,]<-as.numeric(tmp), X[i,]<-as.numeric(tmp[-1]) ) 
      print(i)
   } 
   close(fileIn)
   rownames(X)<-IDx
   colnames(X)<-colNames

   if(weighting){ weight<-scan(weight.file,skip=1) }

   S<-0

 ## Imputing, centering, standarizing and weighting. 
   for(i in 1:ncol(X)){
 	    meanXi <- mean(X[,i],na.rm=TRUE)
        # naive imputation
	    X[,i] <- ifelse(is.na(X[,i]),meanXi, X[,i]) 
	    if(ctr){ X[,i]<-X[,i]-meanXi } #centering
          if(std){ X[,i]<-X[,i]/sd(X[,i]) } # standarizing
          if(weighting){ X[,i]<-X[,i]*weight[i] } # weighting
          S<-S+var(X[,i])
          print(i)
   }

   G<-tcrossprod(X)/S

   if(!is.null(colIDy)){ stopifnot(all(IDs%in%rownames(G))) }

   if(!is.null(colIDy)){
     IDs<-factor(IDs,levels=rownames(G))
     Z<-as.matrix(model.matrix(~IDs-1))
     G<-tcrossprod(tcrossprod(Z,G),Z)
   }

   
   save(G,file='G.rda')

   EVD<-eigen(G)
   rownames(EVD$vectors)<-rownames(G)

   save(EVD,file='EVD.rda')

}








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


