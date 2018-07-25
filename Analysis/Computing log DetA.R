source("R/myRfunctions.R")

#Information about interaction matrix and time series data
mu=seq(-2,2,by=0.2)
sigma=seq(0,1, by =0.1)
d=seq(5, 10, by=0.25)
S=60
N=1000 

###################################################################################
#Computing (using simulations) the log Determinant of A for different mu,sigma,d
###################################################################################


#different combinations of mu, sigma and d
params=expand.grid(mu,sigma,d)
ld=numeric(length(mu)) #empty vector to store logDeterminant

k=length(mu)
seq1=seq(k+1,2*k,by=1)
#creating A from different mu,sigma, d
for(i in seq1)
{
  A=makeASymmetric(mu=params[i,1], sigma=params[i,2], d=params[i,3], S=S)
  ld[i-length(seq1)]=log(det(A))
  print(i)
}

#plotting the logDeterminants of the different A's
plot(x=seq1, ld, col="blue",xlab="combinatioins of mu, sigma,d",ylab="logDetA", 
     title(main="The log Determinant of A for different mu, sigma, d(computed w/simulations"))

avg=mean(ld)#computing the mean of the logDetA

#write.csv(ld, file="testing.tst")



################################################################
#Computing the logDetA (analytically) for different mu,sigma, d
################################################################

# ld2=numeric(nrow(params))#create empty vector
# 
# #solving the analytical solution for log Det A for differerent mu,sigma, d
# for (i in 1:nrow(params))
# {
#   ld2[i]=logDetA(mu=params[i,1], sigma=params[i,2], d=params[i,3], S=S)
#   print(i)
# }
# 
# #plotting logDet A
# plot(x=c(1:nrow(params)), ld2, col="red",xlab="combinatioins of mu, sigma,d",ylab="logDetA", 
#      title(main="The log Determinant of A for different mu, sigma, d(computed analytically)"))
# 
# avg2=mean(ld2)#computing the mean of the logDetA

