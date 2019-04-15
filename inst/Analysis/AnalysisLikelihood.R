##Analysis of Likelihood functions

#A&X_data
source("Data/A&X_data.R")


#myRfunctions
source("R/myRfunctions.R")

Eg=eigen(A)
minEg=min(Re(Eg$values))
maxEg=max(Re(Eg$values))

##Initialization
#
list=seq(-10,10, by=0.1)#values for mu, sigma, d
Y=matrix(rnorm((length(list))*3), nrow=3, ncol=length(list)) #initializing vector for plots, Y

##get from mathematica ##values for plogD w.r.t mu, sigma, d
##
  # Dd=
  # Ddval=
  # Dsigma=
  # Dsigval=sigma
  Dm=1/list
Dmval=1/mu

#Generating Data
A=makeASymmetric(mu=mu, sigma=sigma, d=d, S=S) #generating A
X=makeX(mu=mu, S=S, N=N, A=A) #generating X


#####################################################
##PLOTING PARTIAL DERIVATIVES AGAINST mu, sigma, d
#####################################################





#Y is a list holding the values of partial of logP wrt mu
for(i in 1:length(list))
{
  Y[1,i]=pDm(S, N, X, Dm[i])#for mu=list, sigma=sigma, d=d
  Y[2,i]=pDm(S, N, X, Dmval)#for mu=mu, sigma=list, d=d
  Y[3,i]=pDm(S, N, X, Dmval)#for mu=mu, sigma=sigma, d=list
}

# #Y1 is a list holding the values of partial of logP wrt mu
# for(i in 1:length(list))
# {
#   Y1[1,i]=pDd(N, X, Dd[i])#for mu=mu, sigma=sigma, d=list
#   Y1[2,i]=pDd(N, X, Ddval)#for mu=list, sigma=sigma, d=Ddval
#   Y1[3,i]=pDd(N, X, Ddval)#for mu=mu, sigma=list, d=Ddval
# }

# #Y2 is a list holding the values of partial logP wrt mu
# for(i in 1:length(list))
# {
#   Y2[1,i]=pDsigma(S, N, X, Dsigma[i])#for mu=mu, sigma=list, d=d
#   Y2[2,i]=pDsigma(S, N, X, Dsigval)#for mu=list, sigma=sigma, d=d
#   Y2[3,i]=pDsigma(S, N, X, Dsigval)#for mu=mu, sigma=sigma, d=list
# }


#plotting mu, sigma,d against dLogP/dmu
x=(list) #x values for graph
plot(x, Y[1,], type="p", col="red", xlab="mu=red, sigma,d=green", ylab="dLogP/dmu")
lines(x, y=Y[2, ], type="p", col="blue")
lines(x, y=Y[3, ], type="p", col="green")

# #plotting mu, sigma, d against dLogP/dd
# x=(list) #x values for graph
# plot(x, Y1[1,], type="p", col="red", xlab="d=red, sigma,mu=green", ylab="dLogP/dd")
# lines(x, y=Y1[2, ], type="p", col="blue")
# lines(x, y=Y1[3, ], type="p", col="green")
# 
# #plotting mu, sigma, d against dLogP/dsigma
# x=(list) #x values for graph
# plot(x, Y2[1,], type="p", col="red", xlab="sigma=red, mu,d=green", ylab="dLogP/dmu")
# lines(x, y=Y2[2, ], type="p", col="blue")
# lines(x, y=Y2[3, ], type="p", col="green")



