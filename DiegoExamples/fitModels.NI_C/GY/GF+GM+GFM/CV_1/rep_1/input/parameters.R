folds   <- 1:5 # -999 (no folds,because it is predicting the new combinations of Parents 1126Notes)
 nIter  <- 50000
 burnIn <- 40000
phenotype.file <- '../../../../../../root.Data/Z2_1.csv' 
 AB <- list() 
AB[[1]] <- '../../../../../../G/GcI.P1/output/EVD.rda'       # path to pehnotype file # G1
AB[[2]] <- '../../../../../../G/GcI.P2/output/EVD.rda'       # path to pehnotype file  # G2
AB[[3]] <- '../../../../../../I/GcP1P2/output/EVD.rda'       # path to pehnotype file  # interaction
 type <- c('RKHS','RKHS','RKHS','RKHS','RKHS','RKHS') 
 colENV  <- NULL  # column in phenotype file that gives the id of the environment
# which col for the phenotype environmental
colVAR  <- 3  # column in phenotype file that gives the id of the variety # this is the column for (P1XP2)

 colPhen <- 17
 colCV <- 13


CV0 <- FALSE
 ESC <- FALSE 
 r <- 1
 set.seed(1)
