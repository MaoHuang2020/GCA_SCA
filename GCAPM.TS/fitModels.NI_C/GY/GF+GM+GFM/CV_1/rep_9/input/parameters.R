 folds   <- 1:5
 nIter  <- 50000
 burnIn <- 40000
phenotype.file <- '../../../../../../root.Data/Z2_9.csv' 
 AB <- list() 
 AB[[1]] <- '../../../../../../G/GcI.P1/output/EVD.rda'       # path to pehnotype file 
 AB[[2]] <- '../../../../../../G/GcI.P2/output/EVD.rda'       # path to pehnotype file 
 AB[[3]] <- '../../../../../../I/GcP1P2/output/EVD.rda'       # path to pehnotype file 
 type <- c('RKHS','RKHS','RKHS','RKHS','RKHS','RKHS') 
 colENV  <- NULL  # column in phenotype file that gives the id of the environment 
 colVAR  <- 3  # column in phenotype file that gives the id of the variety 
 colPhen <- 17
 colCV <- 13
 CV0 <- FALSE 
 ESC <- FALSE 
 r <- 1
 set.seed(9)
