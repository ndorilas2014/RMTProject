source("RMTProject/R/myRfunctions.R")

#Information about interaction matrix and time series data
mu=seq(-2,2,by=0.2)
sigma=seq(0,1, by =0.1)
d=seq(5, 10, by=0.25)
S=60
N=1000 

###################################################################################
#Computing (using simulations) the log Determinant of A for different mu,sigma,d
###################################################################################


#different combinations of mu, sigma and d (mu->sigma->d)
params=expand.grid(mu,sigma,d)
#different combinations of sigma, mu,d (sigma->mu->d)
params2=expand.grid(sigma,mu,d)
#different combinations of d, mu, sigma (d->mu->digma)
params3=expand.grid(d,mu,sigma)


seq1=seq(1,length(mu)*length(sigma)*length(d),by=1) #length of ld
k=length(seq1)
ld=numeric(k) #empty vector to store logDeterminant
#creating A from different mu,sigma, d
for(i in seq1)
{
#  A=makeASymmetric(mu=params[i,1], sigma=params[i,2], d=params[i,3], S=S)
#  A=makeASymmetric(mu=params2[i,2], sigma=params2[i,1], d=params2[i,3], S=S)
  A=makeASymmetric(mu=params3[i,2], sigma=params3[i,3], d=params3[i,1], S=S)
  
  ld[i]=(1/S)*log(det(A))
  print(i)
}

#Save plots
jpeg("RMTProject/Data/Plots/The log Determinant of A for different d(simulated).jpeg")
#plotting the logDeterminants of the different A's
plot(x=d, ld[1:length(d)], col="blue",xlab="d->5:10",ylab="logDetA", 
     title(main="The log Determinant of A for different d(computed w/simulations)"))
dev.off()


#Saving values of logDetA to Data folder

# W1=cbind(ld,params)# log determinants and paramater values mu->sigma->d
# write.csv(W1, file="RMTProject/Data/logDeterminants of A/logDetA mu_sigma_d(simulated).txt")

# W2=cbind(ld,params2)# log determinat and parameter values sigma->mu->d
# write.csv(W2, file="RMTProject/Data/logDeterminants of A/logDetA sigma_mu_d(simulated).txt")

# W3=cbind(ld,params3)#log determinant and paramter values d->mu->sigma
# write.csv(W3, file="RMTProject/Data/logDeterminants of A/logDetA d_mu_sigma(simulated).txt")




################################################################
#Computing the logDetA (analytically) for different mu,sigma, d
################################################################
ld2=numeric(k)#create empty vector

#solving the analytical solution for log Det A for differerent mu,sigma, d
for (i in seq1)
{
  ld2[i]=DetA(mu=params[i,1], sigma=params[i,2], d=params[i,3], S=S)
#  ld2[i]=DetA(mu=params2[i,2], sigma=params2[i,1], d=params2[i,3], S=S)
#  ld2[i]=DetA(mu=(params3[i,2]), sigma=(params3[i,3]), d=(params3[i,1]), S=S)
  print(i)
 
}
#Saving log determinants to Data folder

# W1=cbind(ld2,params)# log determinants and paramater values mu->sigma->d
# write.csv(W1, file="RMTProject/Data/logDeterminants of A/logDetA mu_sigma_d(analytical).txt")

# W2=cbind(ld2,params2)# log determinat and parameter values sigma->mu->d
# write.csv(W2, file="RMTProject/Data/logDeterminants of A/logDetA sigma_mu_d(analytical).txt")
# 
# W3=cbind(ld2,params3)#log determinant and paramter values d->mu->sigma
# write.csv(W3, file="RMTProject/Data/logDeterminants of A/logDetA d_mu_sigma(analytical).txt")







#Saving Plot
jpeg("RMTProject/Data/Plots/The log Determinant of A for different sigma(computed analytically).jpeg")
#plotting logDet A
plot(x=sigma, ld2[1:length(sigma)], col="red",xlab="sigma->0:1",ylab="logDetA",
     title(main="The log Determinant of A for different sigma(computed analytically)"))
dev.off()

#################################################################
#read from fies in logdeterminant folder
#################################################################
x1=read.csv("RMTProject/Data/logDeterminants of A/logDetA mu_sigma_d(simulated).txt")
x2=read.csv("RMTProject/Data/logDeterminants of A/logDetA sigma_mu_d(simulated).txt")
x3=read.csv("RMTProject/Data/logDeterminants of A/logDetA d_mu_sigma(simulated).txt")

y1=read.csv("RMTProject/Data/logDeterminants of A/logDetA mu_sigma_d(analytical).txt")
y2=read.csv("RMTProject/Data/logDeterminants of A/logDetA sigma_mu_d(analytical).txt")
y3=read.csv("RMTProject/Data/logDeterminants of A/logDetA d_mu_sigma(analytical).txt")


xval=x1[,2] #x axis (simulated logDeterminants) 
yval=y1[,2] #y axis (analytical logDeterminants) 

# xval=x2[,2]#simulated mu=0, sigma=0.1:1, d=5
# yval=y2[,2]#analytical mu=0, sigma=0.1:1, d=5
# 
# xval=x3[,2]#simulated mu=0, sigma=0.1, d->5:10
# yval=y3[,2]#analytical mu=0, sigma=0.1, d->5:10

plot(xval,yval, xlab="Simulated Value of The logDeterminants", 
     ylab="Analytical Value of The logDeterminants", col="purple",
     title(main="Simulated log Determinants vs Analytical log Determinants"))
#simulated vs. analytical value

#####################################################################
#Simulated logDetA vs Analytically Computed logDetA
#####################################################################

#Saving Plot
jpeg("RMTProject/Data/Plots/Simulated log(DetA) vs. Analytically Computed log(DetA) (d).jpeg")
#ploting logDetA(simulated vs. computed)
plot(ld,ld2, col="purple", xlab="logDetA computed w/simulations", ylab="logDetA computed analytically", 
     title(main="log(DetA) (sigma=0,d->5:10,mu=0.2)"))
dev.off()








