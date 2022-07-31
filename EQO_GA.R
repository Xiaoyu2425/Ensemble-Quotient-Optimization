#' Optimization of microbiome Ensemble Quotient with Genetic Algorithm
#' 
#' This function optimizies Ensemble Quotient of a microbiome 
#' towards trait-based coarse-graining. It is suitable for 
#' a uniform trait, a continuous trait or a categorical trait.
#' 
#' @param pattern A character speficifying the type of trait, including \code{"u"} for a uniform trait, \code{"c"} for a continuous trait and \code{"d"} for a discrete trait
#' @param M A matrix of taxa in the microbiome (samples as rows and taxa as columns)
#' @param y A vector or matrix of trait. y is optional when pattern is \code{"u"}. y is a vector when pattern is \code{"c"}. y is a binary matrix when pattern is \code{"d"} whose rows are samples and columns are categories.
#' @param Nmax A number (default 10) for regularization, specifying the maximal number of taxa allowed in the final group
#' @param amin A number (default 0) specifying the lower bound of average relative abundance of the final group across all samples (only required for a uniform trait to avoid trivial solutions).
#' @param amax A number (default 1) specifying the upper bound of average relative abundance of the final group across all samples (only required for a uniform trait to avoid trivial solutions).
#' @param maxIter A number (default 400) for GA package, specifying the maximum number of iterations to run before the GA search is halted.
#' @param popSize A number (default 100) for GA package, specifying the initial population size for the GA search. 
#' @param parallel Optional argument for GA package (default TRUE), indicating whether parallel computing is needed. See GA documentation for full details.
#' @param monitor Optional argument for GA package (default plot). \code{monitor=plot} indicates that algorithm implementation is visualized as a plot. \code{monitor=gaMonitor} indicates that algorithm implementation is printed as text. \code{monitor=FALSE} indicates the output should be suppressed. 
#' @return A list with the following components.
#' \itemize{
#'   \item x - A binary numeric vector specifying presence/absence of taxa in the final group
#'   \item members - A character vector Names of taxa selected in the final group
#'   \item abundance - A numeric vector with relative abundance of the final group in each sample.
#'   \item y - A number of the group performance. It equals to coefficient of variation when the trait is uniform, for which a smaller value implies higher stability. It eqauls to correlation coefficient when the trait is continuous, for which a larger absolute value implies stronger association. It equals to coefficient of determination when the trait is discrete, for which a larger value implies better discrimination across different categories.
#' }
#' @export

EQ_optim<-function(pattern,M,y,Nmax=10,amin=0,amax=1,popSize=100,maxIter=200,parallel=TRUE,monitor=plot){
	
	if(missing(y)){y<-1;pattern<-"u"}
	
	an<-colMeans(M)
	m<-nrow(M)
	n<-ncol(M)

	c1<-function(x){(rep(1,n) %*% x)-Nmax} # maximal number of taxa
	c2<-function(x){amin-(an %*% x)} # minimal average relative abundance
	c3<-function(x){(an %*% x)-amax} # maximal average relative abundance

	if(pattern=="c"){
		M0<-apply(M,2,function(x){x-mean(x)})
		y0<-y-mean(y)
		P<-t(M0) %*% M0
		Q<-t(M0) %*% y0 %*% t(y0) %*% M0
		Q2<-t(M0) %*% y0
		
		fitness<-function(x){
			pen<-sqrt(.Machine$double.xmax)
			(t(x) %*% Q2)/sqrt((t(x) %*% P) %*% x)-(pen*max(c1(x),0))
		}
	}

	if(pattern=="d"){
		M0<-apply(M,2,function(x){x-mean(x)})
		y0<-y
		L<-matrix(0,nrow=ncol(y0),ncol=ncol(y0))
		diag(L)<-1/sqrt(colSums(y0))
		P<-t(M0) %*% M0
		Q<-t(M0) %*% y0 %*% L %*% L %*% t(y0) %*% M0
		Q2<-t(M0) %*% y0 %*% L
		
		fitness<-function(x){
			pen<-sqrt(.Machine$double.xmax)
			(t(x) %*% Q2)/sqrt((t(x) %*% P) %*% x)-(pen*max(c1(x),0))
		}
	}

	if(pattern=="u"){
		M0<-M
		e<-as.matrix(rep(1,m),ncol=1)
		P<-t(M0) %*% M0 - ((2/n) * (t(M0) %*% e %*% t(e) %*% M0)) + ((1/(n^2)) * (t(M0) %*% e %*% t(e) %*% e %*% t(e) %*% M0))
		Q<-((1/(n^2)) * t(M0) %*% e %*% t(e) %*% M0)
		
		fitness<-function(x){
			pen<-sqrt(.Machine$double.xmax)
			((t(x) %*% Q) %*% x)/((t(x) %*% P) %*% x)-(pen*max(c1(x),0))-(pen*max(c2(x),0))-(pen*max(c3(x),0))
		}
		
	}

	GA<-GA::ga("binary",
		fitness=fitness,
		nBits=n,
		maxiter=maxIter,
		popSize=popSize,
		names=colnames(M0),
		monitor=monitor,
		parallel=parallel)
	
	if(pattern=="c"){
		x<-as.numeric(GA@solution)
		y<-max(GA@fitness[!is.na(GA@fitness)])
		members<-colnames(M0)[GA@solution==1]
		abundance<-rowSums(cbind(rep(0,m),rep(0,m),M[,GA@solution==1]))
		r<-cor(abundance,y0)^2
		return(list(x=x,members=members,abundance=abundance,y=r))
	}

	if(pattern=="d"){
		x<-as.numeric(GA@solution)
		y<-max(GA@fitness[!is.na(GA@fitness)])
		members<-colnames(M)[GA@solution==1]
		abundance<-rowSums(cbind(rep(0,m),rep(0,m),M[,GA@solution==1]))
		s<-abundance-mean(abundance)
		R2<-(t(s) %*% y0 %*% L %*% L %*% t(y0) %*% s)/(t(s) %*% s)
		return(list(x=x,members=members,abundance=abundance,y=R2))
	}

	if(pattern=="u"){
		x<-as.numeric(GA@solution)
		y<-max(GA@fitness[!is.na(GA@fitness)])
		members<-colnames(M0)[GA@solution==1]
		abundance<-rowSums(cbind(rep(0,m),rep(0,m),M0[,GA@solution==1]))
		CV<-sd(abundance)/mean(abundance)
		return(list(x=x,members=members,abundance=abundance,y=CV))
	}	
}
