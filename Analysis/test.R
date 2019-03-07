source("../R/myRfunctions.R")

# 
# 
#  # A=matrix(rnorm(25, mean=0, sd=1), 5, 5)
#  # X=matrix(1,5,5)
#  # sum(A%*%X)
#  # X%*%A%*%t(X)
# DetA(1,5,0.5,20)
# y=c(0.3,7,0.5)
# S=100
# N=100
# A=matrix(rnorm(S*S, mean=y[1], sd=y[3]), S, S)
# ofx(y,S,N)
# 
# 
# #cubintegrate(ofx(y,S=10,N=10),lower=c(-2,5,0), upper=c(2,10,1),fDim=3, method="hcubature")
# 
# 
# 
# 
# muV=seq(-2,2,by=0.2)
# sigmaV=seq(0,1, by =0.1)
# dV=seq(5, 10, by=0.25)
# 
# m=length(muV)*length(sigmaV)*length(dV)
# 
# L=0
# t=numeric(m)
# params=expand.grid(muV,sigmaV,dV)
# 
# seq1=seq(1,m,by=1) #length of ld
# ld=numeric(m) #empty vector to store logDeterminant
# 
# for (i in seq1)
# {
# 
#   A=makeASymmetric(mu=params[i,1], sigma=params[i,2], d=params[i,3],S)
#   X=makeX(mu=params[i,1], S=S, N=N,A=A)
#   
#   
#   ld[i]=sqrt(DetA(params[i,1], params[i,2], params[i,3], S=S))
#   #/((2*pi)^(S/2))*exp((1/2)*sum(X%*%A%*%X))
# 
# }
# 
# 
# 
# 


###Figuring out where error is in P_x

S=20 #matrix size
mu=-1 #mean
sigma=0.5 #deviation
d=10 #diagonal elements
N=20 #number of data points per species
e=-1
muV=seq(-10,10,by=0.5)
sigmaV=seq(0.1,10, by =0.2)
dV=seq(0, 10, by=1)
#ind=numeric(nrow(params))
params=expand.grid(muV,sigmaV,dV)

p=lapply(1:nrow(params), function(i){
    A=makeASymmetric(mu=params[i,1],sigma=params[i,2],d=params[i,3], S=S)
    #X=tryCatch(expr=makeX(mu=params[i,1],sigma=params[i,2],d=params[i,3],S=S,N=N))
    X=try(makeX(mu=params[i,1],sigma=params[i,2],d=params[i,3],S=S,N=N))#createX
    if(class(X)=='try-error'){
        return(NA)
    }else{
    return(P_x(X,A))
    }
    
})

ind=lapply(1:nrow(params),function(i){
    if(is.na(as.numeric((p[i])))){
        print(paste("problem at index: ", i))
        print(params[i,])
        return(i)
    }
    else{
        return(0)
    }
})

# A=makeASymmetric(mu=params[i,1],sigma=params[i,2],d=params[i,3], S=S)
# X=makeX(mu=params[i,1],sigma=params[i,2],d=params[i,3],S,N)
# 
# p[i]=P_x(X,A)

# 
# if(errorND(A)==1)
# {
#     print("Error: The matrix is not diagonally stable")
# } else
#     
# {
#     ###Variables
#     
#     N=500
#     dt=0.01 #time steps
#     y=matrix(0,nrow=S, ncol=N) #population
#     #yint=rnorm(mean=7, sd=3,n=10) #initial condition
#     yint=15
#     
#     
#     
#     y<-createTS(A=A,dt=dt,y0=1,N=N,S=S,sigma=15)
#     
# }

# mu=1
# sigma=0.5
# d=10
# S=10
# N=10
# B=100
# 
# X=makeX(mu,sigma,d,S,N)
# optim(par=c(mu,sigma,d),fn=MonteCarlo,X=X,B=B,method="BFGS")
