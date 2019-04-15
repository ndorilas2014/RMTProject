################
#RMT Results
################

source("RMTProject/R/myRfunctions.R")

mu=0
sigma=1
rhoB=0.5
d=-5
S=500

A=matrix(rnorm(S*S, mean=mu, sd=sigma), S, S)#creating random matrix, random entries

#fill off diagonal entries of A with random samples from bivariate dis. of pairs(Mij, Mji)
for (i in (1:(S-1)))
{
  for (j in (1+i):(S))
  {
    B=rnorm(1, (mu), (sigma))#normal distribution, mean=mu, deviation=sigma
    A[i,j]=B[1] #creating a symmetric matrix by giving aij and aji the same value
    A[j,i]=B[1] #^^^^^^
  }
}

##bivariate distribution
A=matrix(rnorm(S), S, S)
Sig=matrix(c(sigma^2, rhoB*sigma^2, rhoB*sigma^2, sigma^2), 2) #covariance matrix

#fill off diagonal entries of A with random samples from bivariate dis. of pairs(Mij, Mji)
for (i in (1:(S-1)))
{
  for (j in (1+i):(S))
  {
    B=rnorm(1, mu, sigma) #bivariate distribution, 1 sample
    A[i,j]=B[1]
    A[j,i]=B[1]
  }
}

for (i in 1:S)
{
  A[i,i]=-d
}


eg=eigen(A)
yval=eg$values

plot(yval, xlab="real axis", ylab="imaginary axis")
