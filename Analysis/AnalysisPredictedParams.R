###########################################
#PLOTTING ACTUAL VALUES OF mu, sigma, d 
#AGAINST PREDITED VALUES of mu,sigma d
###########################################


#Functions
source("RMTProject/R/myRfunctions.R")


#Information about interaction matrix and time series data
mu=seq(-4,7,by=0.25)
sigma=seq(0,1, by =0.1)
d=seq(5, 10, by=0.25)
S=60
N=100 

#different combinations of mu, sigma and d
params=expand.grid(mu,sigma,d)

#to store the predicted values of mu,sigma,d


q=matrix(0, N, 3)#holds predicted values of mu,sigma, d
