source("RMTProject/R/myRfunctions.R")

#Information about interaction matrix and time series data
mu=seq(0.2,2,by=0.2)
sigma=seq(0.1,1, by =0.1)
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

k=length(mu)
seq1=seq(1,length(mu),by=1) #length of ld
ld=numeric(length(seq1)) #empty vector to store logDeterminant
#creating A from different mu,sigma, d
for(i in seq1)
{
  A=makeASymmetric(mu=params[i,1], sigma=params[i,2], d=params[i,3], S=S)
#  A=makeASymmetric(mu=params2[i,2], sigma=params2[i,1], d=params2[i,3], S=S)
#  A=makeASymmetric(mu=params3[i,2], sigma=params3[i,3], d=params3[i,1], S=S)
  
  ld[i]=log(det(A))
  print(i)
}

#Save plots
jpeg("RMTProject/Data/Plots/The log Determinant of A for different d(simulated).jpeg")
#plotting the logDeterminants of the different A's
plot(x=mu, ld, col="blue",xlab="mu->0.2:2",ylab="logDetA", 
     title(main="The log Determinant of A for different mu(computed w/simulations)"))
dev.off()


#Saving values of logDetA to Data folder

#write.csv(ld, file="logDetA for varying mu, fixed sigma,d(simulated).txt")
#write.csv(ld, file="logDetA for varying sigma, fixed mu,d(simulated).txt")
#write.csv(ld, file="logDetA for varying d, fixed mu,sigma(simulated).txt")




################################################################
#Computing the logDetA (analytically) for different mu,sigma, d
################################################################

ld2=numeric(length(seq1))#create empty vector

#solving the analytical solution for log Det A for differerent mu,sigma, d
for (i in seq1)
{
  ld2[i]=DetA(mu=(params[i,1]), sigma=(params[i,2]), d=(params[i,3]), S=S)
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
     title(main="Simulated log(DetA) vs. Analytically Computed log(DetA) (mu->0.2:2)"))
dev.off()

