##libraries
library(matlib)
library(MASS)
library(mvtnorm)
library(rootSolve)
library(cubature)
library(Matrix)
library(Brobdingnag)


##Functions
#####################################################################################




#population model
TS<-function(A, dt, y, xi, i,j)
{
  y[i,j]=y[i,(j-1)]+dt*((A%*%y[,(j-1)])[i]+xi)
  return(y)
}

createTS<-function(A,dt,y0,N,S,sigma)
{
  y=matrix(0,nrow=S, ncol=N) #population
  y[,1]=y0 #initial condition
  x1=rnorm(n=S, mean=0, sd=sigma)
  x2=matrix(rnorm(n=N*S, mean=0, sd=sigma),S,N)
  
  
  for(i in 1:(N-1))
  {
    
    y[,i+1]=y[,i]+dt*(A%*%y[,i]+x2[,i+1])
  
  }
  return(y)
  
}





########################################################################################

##Partial Derivative Functions

#Returns the partial derivative of the rhs logP w.r.t mu 
rhsmu<-function(S, N, X)
{
  sy=numeric(N) #creates empty vector of size N to store the squared sums of columns
  #stores the sum of all columns of X, then squares it
  for (a in 1:N)
  {
    sy[a]=(sum(X[,a]))^2
  }
  SY=sum(sy) #sum of the sums
  
  return ((SY/(S*N))) #pLogP wrt mu
}

#partial derivative of logDet wrt mu
pDm<-function(mu,d)
{
  return(1/(d+mu))
}

#Returns the rhs of partial derivative of logP w.r.t d
rhsD<-function(S, N, X)
{
  sy=numeric(N) #creates empty vector of size N
  #stores the sum of all columns of X, then squares it
  for (a in 1:N)
  {
    sy[a]=sum((X[,a])^2)
  }
  SY=sum(sy) # sum of the sums
  
  sy2=numeric(N) #creates empty vector of size N
  #stores the sum of all columns of X, then squares it
  for (a in 1:N)
  {
    sy2[a]=(sum(X[,a]))^2
  }
  SY2=sum(sy2) # sum of the sums
  
  return ((-SY2/S^2/N)+(SY/(S*N))) 
}

#partial of DetA wrt D
pDetD<-function(sigma,d)
{
  x=sigma^2/d^2
  f=1/(1+sqrt(1-4*x)-2*x)
  
  return((1/d)+(2*sigma^2/d^3)*f)
}

#Returns the rhs of partial derivative of logP w.r.t sigma
rhssigma<-function(S, N, X)
{
  # ab=expand.grid(1:N, 1:N) #creates data frame#one row for each combination of a,b, 1-N
  # foo=numeric(nrow(ab)) # creates empty vector of length number of rows in ab
  # #foo stores the matrix multiplication of sum_i(sum_ab(X[,a])(X[,b]))
  # for(i in (1:nrow(ab)))
  # {
  #   foo[i]=((X[,ab[i,1]])%*%(X[,ab[i,2]]))^2 #multiplies all combinations of x[,a] and x[,b]
  # }
  # SY=sum(foo)#stores the sum of this new vector
  # 
  # return ((SY/(S*N*2)/S/1+N/S)) #rhs pLogp wrt sigma
  return(sum((X %*% t(X))^2 ) /(2*S*N) / S /(1+N/S))
}

#partial detA wrt sigma
pDetS<-function(sigma,d)
{
  b=sigma^2/d^2
  return(1./(d^2) /( 1+ sqrt(1-4*b)-2*b ))
}




#################################################################################

##Optimized mu,sigma,d

# #Returns optimized mu
# MuF<-function(S,N,X)
# {
#   sy=numeric(N) #creates empty vector of size N to store the squared sums of columns
#   #stores the sum of all columns of X, then squares it
#   for (a in 1:N)
#   {
#     sy[a]=(sum(X[,a]))^2
#   }
#   SY=sum(sy) #sum of the sums
#   
#   return(-(S*N)/(SY)) #mu
# }

#HypergeometricPFQ function
hypergeo3F2 <- function(x){
  
  k <- 0
  s <- 0
  oldldx <- 0
  ldx <- 1
  dx <- 1
  if(x==0){
    return(1)
  }else{
  while( dx / s > 0.0001 ){
    ldx <- lgamma(1.5+k) + log(4) - log(1+k) - 0.5*log(pi) -
      lgamma(3+k) + k * log(x)
    dx <- exp(ldx)
    s <- s + dx
    oldldx <- ldx
    k <- k + 1
  }
  }
  
  return(s)
  
}

#NEEDS WORK
#return functions for optimizing d and sigma
Opfunct<-function(x, S, N, X) {
  
  x[1]=exp(x[1])
  x[2]=exp(x[2])
  #X=makeX(mu=x[3], S=S, N=N, A=makeASymmetric(mu=x[3],sigma=x[1],d=x[2], S=S))
  return(c( F1=-rhssigma(S,N,X)+pDetS(sigma=x[1],d=x[2]),
            F2=-rhsD(S,N,X)+pDetD(sigma=x[1],d=x[2]),
            F3=-rhsmu(S,N,X)+pDm(mu=x[3],d=x[2])))
}



model <- function(x,S,N,X) c(F1=x[1]+x[2]-13, 
                         F2=2*x[1]-3*x[2]+5)


ofx<-function(x, S, N){
  mu1=0
  mu<-x[1]
  sigma<-x[3]
  d<-x[2]
  B=makeASymmetric(mu=mu1,sigma=sigma,d=d, S=S)
  X=makeX(mu=mu, S=S, N=N, A=B)
  s=numeric(S)
  
  for (i in 1:S) 
  {
    s[i]=sum(X[,i])^2
  }
  
  D=exp(DetA(mu=mu, d=d, sigma=sigma, S=S))
  P1=(-d/2)*sum(X^2)
  P2=sum(s)*(mu/(2))
  P3=(1/2)*sum(X%*%B%*%t(X)) 
  P4=(1/(2*sigma^2))*sum((B^2))
  OF<-exp(P1+P2+P3-P4)*sqrt(1/(2*pi))*(1/sigma)
  return(OF)
}


################################################################################

##TESTING


#function for computing the logDeterminant of A
DetA<-function(mu,d,sigma,S)
{
  I=pi/2 - (sigma/8)*((hypergeo3F2((4*sigma^2)/d^2) - 8*log(d/(2*sigma))))
  I=log(d)-(sigma^2)/(2*d^2)*(hypergeo3F2((4*sigma^2)/d^2))
  return((1/S)*log(d+mu)+((S-1)/S)*I)
  
}






################################################################################

##Functions for Generating data



##Generating species data from P(x|A)=sqrt(detA)/((2*pi)^S/2) * exp((1/2)sum(xiAijxj))
##
makeX<-function(mu, sigma, d, S, N)
{
  # browser()
  mulist=rep(mu, S)#Generating means for X
  A=makeASymmetric(mu,sigma,d,S)
  
  X=t(rmvnorm(N, (mulist/S), inv(A), method="svd")) #Generates X/data from a mv distribution
  #N=ncol(X) #for when X is given not created from mv distribution
  
  return(X)
}



##Generating A from Bivariate Distribution
makeA_BV<-function(mu, sigma, rhoB, d, S)
{
  
  ##bivariate distribution
  A=matrix(rnorm(S), S, S)
  Sig=matrix(c(sigma^2, rhoB*sigma^2, rhoB*sigma^2, sigma^2), 2) #covariance matrix
  
  #fill off diagonal entries of A with random samples from bivariate dis. of pairs(Mij, Mji)
  for (i in (1:(S-1)))
  {
    for (j in (1+i):(S))
    {
      B=mvrnorm(1, c(mu,mu), Sig) #bivariate distribution, 1 sample
      A[i,j]=B[1]
      A[j,i]=B[2]
    }
  }
  
  #fill diagonal entries of A with -d
  for (i in 1:S)
  {
    A[i,i]=-d
  }
  
  A=A/sqrt(S) #scale matrix by sqrt(S)
  
  return(A)
}



##Generating Symmetric A
#Information about A-interaction matrix and X-time series data--mu, sigma, d, S

makeASymmetric<-function(mu, sigma, d, S)
{
  ##GENERATING A

  A=matrix(rnorm(S*S, mean=mu, sd=sigma), S, S)#creating random matrix, random entries
  R=rnorm(((S*S)-S)/2, (mu/S),(sigma/sqrt(S)))

  A=forceSymmetric(A)

  diag(A)=d


  return(as.matrix(A))
}

P_x<-function(X, A)
{
    # note: we're always returning the log of P(X|A)
    
    S=nrow(A)
    N=ncol(X)
    
    logDA = determinant(A, logarithm = TRUE)
    logDA = logDA$modulus * logDA$sign
    
    out <- N/2 * (logDA - S * log(2 * pi)) - 1/2 * sum(diag(t(X)%*%A%*%X))
    
    return(as.numeric(out))
}

MonteCarlo<-function(x,X, B)
{
    mu=x[1]
    sigma=x[2]
    d=x[3]
    
    sumB=0
    N=ncol(X)
    S=nrow(X)
    

    L=lapply(1:B, function(i){
        A=makeASymmetric(mu, sigma, d, S=S)
        return(P_x(X, A))
    })
    L=mean(unlist(L),na.rm=TRUE)
    #browser()
    
    return(L)
    
}

RMToptimize<-function(X,niter=100)
{
    mu=1
    sigma=0.5
    d=10
    S=nrow(X)
    N=ncol(X)
    B=niter
    
    #X=makeX(mu,sigma,d,S,N)
    out<-optim(par=c(mu,sigma,d),fn=MonteCarlo,X=X,B=B, method = 'L-BFGS-B', 
          lower = c(-10, 0, 0), upper = c(10, 10, 100),
          control = list(fnscale = -1))
    return(out)
    
}


##Error if A is not positive definite
errorPD<-function(A)
{
  Eg=eigen(A)
  minEg=min(Re(Eg$values))
  
  if(minEg<=0)
  {
    return(1)
  }else{
    return(0)
  }
}





########################################################################

##Error Messsages



#Error is A is not negative definite
errorND<-function(A)
{
  Eg=eigen(A)
  maxEg=max(Re(Eg$values))
  
  if(maxEg>=0)
  {
    return(1)
  }else{
    return(0)
  }  
}