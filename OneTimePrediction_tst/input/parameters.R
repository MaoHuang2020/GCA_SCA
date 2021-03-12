###Paremeters

phenotype.file <- '/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction/data/Y.csv'          # path to pehnotype file
mm.file        <- '/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction/data/A.csv'          # path to covariates file
weight.file     <- NULL    # path to weight file if NULL no weights are used


ctr 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
std 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
weighting 	<- !is.null(weight.file) 	# Valid only if mm.file is not NULL


colIDy <- 3  # column in phenotype file that gives the IDs that link observations to covariates or grouping factor
# or NULL if ID is not used
# 2 #(basicallythe P1 column ---- 1126Notes)


# folds   <- -999 # -999 (no folds,because it is predicting the new combinations of Parents 1126Notes)
# nIter  <- 5000
# burnIn <- 4000
# phenotype.file <- '/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction/data/CrossMergeDwPMTrial1.csv' 
# AB <- list() 
# AB[[1]] <- '/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction/output/EVD.rda'       # path to pehnotype file # G1
# AB[[2]] <- '/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction/output/EVD.rda'       # path to pehnotype file  # G2
# AB[[3]] <- '/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction/output/EVD.rda'       # path to pehnotype file  # interaction
# type <- 'RKHS'   #'RKHS','RKHS','RKHS','RKHS','RKHS' 
# colENV  <- 6     # column in phenotype file that gives the id of the environment
# # which col for the phenotype environmental
# colVAR  <- 2  # column in phenotype file that gives the id of the variety # this is the column for (P1XP2)
# 
# colPhen <- 5
# #colCV <- 13
# 
# 
# CV0 <- FALSE
# ESC <- FALSE 
# r <- 1
# set.seed(1)
