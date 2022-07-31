# Ensemble Quotient optimization by Boolean Least Square 
# (specifically for a continuous phenotypic variable)

# M, A matrix of taxa in the microbiome (samples as rows and taxa as columns)
# y, A vector for a continuous phenotypic variable
# maxN, penalty term specifying the maximal number of taxa allowed in the selected group
# ub, a sufficiently large number for linearization transformation
# adjust, a sufficiently large number to guarantee numerically postive semi-definite

EQ_BLS<-function(M,y,maxN,ub,adjust){
	n<-ncol(M)
	M_aug<-cbind(1,M)
	
	L= t(M_aug)%*%M_aug
	L=L/adjust 
	b= -2*(t(M_aug)%*%y)

	model          <- list()
	model$A        <- rbind(c(rep(0,2*n+1), 1),c(0,rep(0,n),rep(1,n),0))
	model$rhs      <- c(ub,maxN) # upper bound on u
	model$sense    <- c('<','<')
	model$obj      <- (1/adjust)*c(b, rep(0, n+1))
	model$Q        <- matrix(0,nrow=2*n+2, ncol=2*n+2)
	model$Q[1:(n+1),1:(n+1)]<-L
	model$vtype    <- c(rep('C', n+1), rep('B', n), 'C')

	qcs = list()
	for (i in 1:n) {
  		qcs[[i]]       <- list()
  		qcs[[i]]$Qc    <- Matrix::spMatrix(2*n+2, 2*n+2, c(n+i+1), c(2*n+2), c(-1.0))
  		qcs[[i]]$rhs   <- 0
  		qcs[[i]]$q     <- c(rep(0, i), 1, rep(0, 2*n+1-i))
  		qcs[[i]]$sense <- '='
	}	

	model$quadcon <- qcs
	result <- gurobi::gurobi(model)

	members<-which(result$x[-1][1:n]>0)
	
	who<-colnames(M)[members]
	abundance<-rowSums(cbind(rep(0,nrow(M)),rep(0,nrow(M)),M[,members]))
	r<-cor(abundance,ext)
	
	return(list(r=r,who=who,abundance=abundance))
}
