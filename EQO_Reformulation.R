# Ensemble Quotient optimization by reformulating into 
# an equivalent mixed integer linear programming problem (Gaur and Arora, 2008)

# M, A matrix of taxa in the microbiome (samples as rows and taxa as columns)
# y, A vector for a continuous phenotypic variable
# ub, a sufficiently large number for linearization transformation
# cutoff, minimal average relative abundance of the group to  avoid generating an empty group or a group with all taxa that are numerically stable but ecologically trivial (only required for a uniform phenotypic variable) 

M<-Microbiome
y<-trait
ub<-100
n<-ncol(M)
cutoff<-0

diagM<-function(x,N){
	dm<-matrix(0,nrow=N,ncol=N)
	diag(dm)<-x
	return(dm)
}

# uniform variable
M0<-M
e<-as.matrix(rep(1,nrow(M0)),ncol=1)
P<-t(M0) %*% M0 - ((2/n) * (t(M0) %*% e %*% t(e) %*% M0)) + ((1/(n^2)) * (t(M0) %*% e %*% t(e) %*% e %*% t(e) %*% M0))
Q<-((1/(n^2)) * t(M0) %*% e %*% t(e) %*% M0)

# continuous variable 
M0<-apply(M,2,function(x){x-mean(x)})
y0<-y-mean(y)
P<-t(M0) %*% M0
Q<-t(M0) %*% y0 %*% t(y0) %*% M0

# categorical variable
M0<-apply(M,2,function(x){x-mean(x)})
y0<-y
L<-matrix(0,nrow=ncol(y0),ncol=ncol(y0))
diag(L)<-1/sqrt(colSums(y0))
P<-t(M0) %*% M0
Q<-t(M0) %*% y0 %*% L %*% L %*% t(y0) %*% M0

# solving the optimization with Gurobi optimizer 
m1<-max(rowSums(abs(P)))
m2<-max(rowSums(abs(Q)))

A<-matrix(0,nrow=(9*n)+4,ncol=(6*n)+1)

A[1:n,1:n] = diagM(1,n)
A[1:n,(n+1):(2*n)] = diagM(1,n)
A[1:n,((4*n)+1):(5*n)] = -P
A[1:n,(6*n)+1] = -m1

A[(n+1):(2*n),1:n] = diagM(1,n)
A[(n+1):(2*n),((4*n)+1):(5*n)] = diagM(2*m1,n)
A[(n+1):(2*n),(6*n)+1] = (-2)*m1

A[((2*n)+1):(3*n),(n+1):(2*n)] = diagM(1,n)
A[((2*n)+1):(3*n),((4*n)+1):(5*n)] = diagM((-2)*m1,n)

A[((3*n)+1):(4*n),((2*n)+1):(3*n)] = diagM(1,n)
A[((3*n)+1):(4*n),((3*n)+1):(4*n)] = diagM(1,n)
A[((3*n)+1):(4*n),((4*n)+1):(5*n)] = -Q
A[((3*n)+1):(4*n),(6*n)+1] = -m2

A[((4*n)+1):(5*n),((2*n)+1):(3*n)] = diagM(1,n)
A[((4*n)+1):(5*n),((4*n)+1):(5*n)] = diagM(2*m2,n)
A[((4*n)+1):(5*n),(6*n)+1] = (-2)*m2

A[((5*n)+1):(6*n),((3*n)+1):(4*n)] = diagM(1,n)
A[((5*n)+1):(6*n),((4*n)+1):(5*n)] = diagM((-2)*m2,n)

A[(6*n)+1,((3*n)+1):(4*n)] = 1
A[(6*n)+1,((4*n)+1):(5*n)] = -m2

A[((6*n)+2):((7*n)+1),((4*n)+1):(5*n)] = diagM(1,n)
A[((6*n)+2):((7*n)+1),((5*n)+1):(6*n)] = diagM(-ub,n)
A[((6*n)+2):((7*n)+1),(6*n)+1] = -1

A[((7*n)+2):((8*n)+1),((4*n)+1):(5*n)] = diagM(1,n)
A[((7*n)+2):((8*n)+1),(6*n)+1] = -1

A[((8*n)+2):((9*n)+1),((4*n)+1):(5*n)] = diagM(1,n)
A[((8*n)+2):((9*n)+1),((5*n)+1):(6*n)] = diagM(-ub,n)

A[((9*n)+2),((6*n)+1)]<-1

A[((9*n)+3),((5*n)+1):(6*n)]<-colMeans(M0)
A[((9*n)+4),((5*n)+1):(6*n)]<-colMeans(M0)

L<-rep(0,(6*n)+1)
L[(n+1):(2*n)] = 1
L[((4*n)+1):(5*n)] = -m1 

model<-list()
model$obj<-L
model$A<-A
model$rhs<-rep(0,nrow(A))
model$rhs[(6*n)+1]<-1
model$rhs[((6*n)+2):((7*n)+1)]<-(-ub)
model$rhs[((9*n)+2)]<-ub
model$rhs[((9*n)+3)]<-cutoff
model$rhs[((9*n)+4)]<-(1-cutoff)
model$sense<-c(rep("=",n),rep("<",n),rep("<",n),rep("=",n),rep("<",n),rep("<",n),"=",rep(">",n),rep("<",n),rep("<",n),"<",">","<")
model$lb<-0
model$vtype<-rep("C",ncol(A))
model$vtype[((5*n)+1):(6*n)]<-"B"

out<-gurobi::gurobi(model)
out$x[(5*n+1):(length(out$x)-1)] # final group selected


