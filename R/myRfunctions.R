##libraries
library(matlib)
library(MASS)
library(mvtnorm)


##Functions
#####################################################################################




#population model
TS<-function(A, dt, y, xi, i,j)
{
  y[i,j]=y[i,(j-1)]+dt*((A%*%y[,(j-1)])[i]+xi)
  return(y)
}





########################################################################################

##Partial Derivative Functions

#Returns the partial derivative of logP w.r.t mu 
pDm<-function(S, N, X, Dm)
{
  sy=numeric(N) #creates empty vector of size N to store the squared sums of columns
  #stores the sum of all columns of X, then squares it
  for (a in 1:N)
  {
    sy[a]=(sum(X[,a]))^2
  }
  SY=sum(sy) #sum of the sums
  
  return (Dm+(SY/(S*N))) #pLogP wrt mu
}

#Returns the parial derivative of logP w.r.t d
pDd<-function(N, X, Dd)
{
  sy=numeric(N) #creates empty vector of size N
  #stores the sum of all columns of X, then squares it
  for (a in 1:N)
  {
    sy[a]=(sum((X[,a])^2))
  }
  SY=sum(sy) # sum of the sums
  
  return (Dd-(SY/N)) #pLogP wrt d
}

#Returns the partial derivative of logP w.r.t sigma
pDsigma<-function(S, N, X, Dsigma)
{
  ab=expand.grid(1:N, 1:N) #creates data frame#one row for each combination of a,b, 1-N
  foo=numeric(nrow(ab)) # creates empty vector of length number of rows in ab
  #foo stores the matrix multiplication of sum_i(sum_ab(X[,a])(X[,b]))
  for(i in (1:nrow(ab)))
  {
    foo[i]=((X[,ab[i,1]])%*%(X[,ab[i,2]]))^2 #multiplies all combinations of x[,a] and x[,b]
  }
  SY=sum(foo)#stores the sum of this new vector
  
  return (Dsigma+(SY/(S*N*2))) #pLogp wrt sigma
}






#################################################################################

##Optimized mu,sigma,d

#Returns optimized mu
MuF<-function(S,N,X)
{
  sy=numeric(N) #creates empty vector of size N to store the squared sums of columns
  #stores the sum of all columns of X, then squares it
  for (a in 1:N)
  {
    sy[a]=(sum(X[,a]))^2
  }
  SY=sum(sy) #sum of the sums
  
  return(-(S*N)/(SY)) #mu
}

#HELP
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
Opfunct<-function(x, params)
{

  #sigma=x[1]
  #d=x[2]

  S=params[1]
  N=params[2]
  X=params[3]

  #H=hypergeo3F2((4*x[1]^2)/(x[2]^2)) #HEEEEEEEEELLLPP
  H=1

  #creating RHS of dlogP wrt d
  sy1=numeric(N) #creates empty vector of size N
  #stores the sum of all columns of X, then squares it
  for (a in 1:N)
  {
    sy1[a]=(sum((X[,a])^2))
  }
  SY1=sum(sy1) # sum of the sums

  #derivative from mathematica
  LHSd=(1/8)*pi*x[1]*((8/x[2])-8*(x[1]^2)*((x[2]^2)*(-(x[2]^2)+2*(x[1]^2)+(x[2]^2)*sqrt(((x[2]^2)
  -4*(x[1]^2))/(x[2]^2)))/(2*(x[1]^4))+H)/(x[2]^3)+(8*H/(x[2]^3)))

  RHSd=(-1)*(SY1/N)

  #creating RHS for dlogP wrt sigma
  ab=expand.grid(1:N, 1:N) #creates data frame#one row for each combination of a,b, 1-N
  foo=numeric(nrow(ab)) # creates empty vector of length number of rows in ab
  #foo stores the matrix multiplication of sum_i(sum_ab(X[,a])(X[,b]))
  for(i in (1:nrow(ab)))
  {
    foo[i]=((X[,ab[i,1]])%*%(X[,ab[i,2]]))^2 #multiplies all combinations of x[,a] and x[,b]
  }
  SY2=sum(foo)#stores the sum of this new vector

  #derivative from mathematica
  LHSsigma= (1/8)*pi*x[1]*(-(8/x[1])+8*x[1]*((x[2]^2)*(-(x[2]^2)+2*(x[1]^2)+(x[2]^2)*sqrt(((x[2]^2)
  -4*(x[1]^2))/(d^2)))/(2*(x[1]^4))+H)/(d^2)-(8*x[1]*H/(d^2)))
  +(1/8)*pi*((-4*(x[1]^2)*H)/(x[2]^2)+8*log(x[2]/(2*x[1])))

  RHSsigma=SY2/(S*N)

  f[1]=LHSsigma+RHSsigma
  f[2]=LHSd+RHSd

  return(f)

}






################################################################################

##TESTING

#NEEDS WORK
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
makeX<-function(mu, S, N, A)
{
  
  mulist=rep(mu, S)#Generating means for X
  
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
  ##
  A=matrix(rnorm(S*S, mean=mu, sd=sigma), S, S)#creating random matrix, random entries
  
  #fill off diagonal entries of A with random samples from bivariate dis. of pairs(Mij, Mji)
  for (i in (1:(S-1)))
  {
    for (j in (1+i):(S))
    {
      B=rnorm(1, (mu/S), (sigma/sqrt(S))) #normal distribution, mean=mu, deviation=sigma
      A[i,j]=B[1] #creating a symmetric matrix by giving aij and aji the same value
      A[j,i]=B[1] #^^^^^^
    }
  }
  
  #setting off diagonal elements to -d
  for (i in 1:S)
  {
    A[i,i]=d 
  }
  
  #A=A/sqrt(S)#resizing 
  
  return(A)
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
    print("A is not negative definite")
    return(1)
  }else{
    return(0)
  }  
}