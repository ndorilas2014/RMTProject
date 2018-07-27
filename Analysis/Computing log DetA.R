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

k=21
seq1=seq(1,length(mu)*length(sigma)*length(d),by=1) #length of ld
ld=numeric(k) #empty vector to store logDeterminant
#creating A from different mu,sigma, d
for(i in seq1[232:252])
{
#  A=makeASymmetric(mu=params[i,1], sigma=params[i,2], d=params[i,3], S=S)
#  A=makeASymmetric(mu=params2[i,2], sigma=params2[i,1], d=params2[i,3], S=S)
  A=makeASymmetric(mu=params3[i,2], sigma=params3[i,3], d=params3[i,1], S=S)
  
  ld[i-231]=log(det(A))
  print(i)
}

#Save plots
jpeg("RMTProject/Data/Plots/The log Determinant of A for different d(simulated).jpeg")
#plotting the logDeterminants of the different A's
plot(x=mu, ld, col="blue",xlab="mu->0.2:2",ylab="logDetA", 
     title(main="The log Determinant of A for different mu(computed w/simulations)"))
dev.off()


#Saving values of logDetA to Data folder

#W1=cbind(ld,params)# log determinants and paramater values mu->sigma->d
#write.csv(W1, file="RMTProject/Data/logDeterminants of A/logDetA mu_sigma_d(simulated).txt")

#W2=cbind(ld,params2)# log determinat and parameter values sigma->mu->d
#write.csv(W2, file="RMTProject/Data/logDeterminants of A/logDetA sigma_mu_d(simulated).txt")

#W3=cbind(ld,params3)#log determinant and paramter values d->mu->sigma
#write.csv(W3, file="RMTProject/Data/logDeterminants of A/logDetA d_mu_sigma(simulated).txt")




################################################################
#Computing the logDetA (analytically) for different mu,sigma, d
################################################################
k=21
ld2=numeric(k)#create empty vector

#solving the analytical solution for log Det A for differerent mu,sigma, d
for (i in seq1[232:252])
{

  ld2[i-231]=DetA(mu=(params3[i,2]), sigma=(params3[i,3]), d=(params3[i,1]), S=S)
  print(i)
 
}
#Saving Plot
jpeg("RMTProject/Data/Plots/The log Determinant of A for different d(computed analytically).jpeg")
#plotting logDet A
plot(x=d, ld2, col="red",xlab="mu->0.2:2",ylab="logDetA",
     title(main="The log Determinant of A for different d(computed analytically)"))
dev.off()






#####################################################################
#Simulated logDetA vs Analytically Computed logDetA
#####################################################################

#Saving Plot
jpeg("RMTProject/Data/Plots/Simulated log(DetA) vs. Analytically Computed log(DetA) (d).jpeg")
#ploting logDetA(simulated vs. computed)
plot(ld,ld2, col="purple", xlab="logDetA computed w/simulations", ylab="logDetA computed analytically", 
     title(main="log(DetA) (sigma=0,d->5:10,mu=0.2)"))
dev.off()








