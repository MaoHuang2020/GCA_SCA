#############################################
 rm(list=ls())

## The working directory should be the analysis directory
 setwd('../output')

## Parameters
 source('../input/parameters.R')

## Checks parameters
 cat('=1===> Reading and Checking Parameters');cat('\n')
 source('../input/parameters.R') 
   if(is.null(G1.file) | is.null(G2.file)){ cat('Please provide path and filename G matrices 1 and 2') }
   cat('The Eigen value decomposition will be perfomed with components \n')

load(G1.file)
G1 <- G
load(G2.file)
G2 <- G

GI <- G1*G2
EVD <- eigen(GI)

save(GI,file='../GP1P2/G.rda')
save(EVD,file='../GP1P2/EVD.rda')

# This is plotting script
pdf(file='Eigen.pdf')
  plot(EVD$vectors[,1],EVD$vectors[,2],main='',xlab='First Component',ylab='Second Component')
  plot(EVD$vectors[,1],EVD$vectors[,3],main='',xlab='First Component',ylab='Third Component')
  plot(EVD$vectors[,2],EVD$vectors[,3],main='',xlab='Second Component',ylab='Third Component')
  plot(EVD$values, main='',xlab='Components',ylab='Eigen-value')
  plot(cumsum(EVD$values)[1:I(min(dim(GI)[1],100))]/sum(EVD$values),main='',xlab='components',ylab='Cumulative Variance')
  abline(h=.8,lty=3)
  abline(v=which(  floor(cumsum(EVD$values)[1:100]/sum(EVD$values)*10) ==8)[1], lty=3)
  dev.off()

quit(save='no')


