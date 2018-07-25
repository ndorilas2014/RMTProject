###########################################
#PLOTTING ACTUAL VALUES OF mu, sigma, d 
#AGAINST PREDITED VALUES of mu,sigma d
###########################################


#Functions
source("R/myRfunctions.R")


#Information about interaction matrix and time series data
mu=seq(-4,7,by=0.25)
sigma=seq(0,1, by =0.1)
d=seq(5, 10, by=0.25)
S=60
N=1000 

#different combinations of mu, sigma and d
params=expand.grid(mu,sigma,d)

#to store the predicted values of mu,sigma,d

Mu1=numeric(length(mu)) #to store predicted value of mu
Sigma1=numeric(length(sigma)) #to store predicted value of sigma
D1=numeric(length(d))# to store predicted value of D1

#stores predicted values for Mu for constant sigma and d
for(r in 1:length(mu))
{
  #create A and X
  
  A=makeASymmetric(mu=params[r,1], sigma=params[r,2], d=params[r,3], S=S) #makeSymmetricA
  #based on params
  
  X=makeX(mu=params[r,1], S=S, N=N, A=A) #generate specie data based on params

  
  #Check if A is positive definite 
  if(errorPD(A)==1)
  {
    print("Matrix is not diagonally stable")
  }else{
    Mu1[r]=MuF(S=S, N=N, X=X) #optimized value of Mu
  }
  
  print(r)


#plotting predicted value against actual value
plot(mu, Mu1, xlab="actual value of Mu", ylab="predicted value of Mu", col="red")

}




#plot(var2, esigma, xlab="actual value of sigma", ylab="predicted value of sigma", col="blue")
#plot(var3, ed, xlab="actual value of d", ylab="predicted value of d", col=green)