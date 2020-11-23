

phenotype.file <- '../../../root.Data/Z2_1.csv'          # path to pehnotype file
mm.file        <- NULL          # path to covariates file
weight.file     <- NULL    # path to weight file if NULL no weights are used


ctr 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
std 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
weighting 	<- !is.null(weight.file) 	# Valid only if mm.file is not NULL
 

colIDy <- 2  # column in phenotype file that gives the IDs that link observations to covariates or grouping factor
             # or NULL if ID is not used
