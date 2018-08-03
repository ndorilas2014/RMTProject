source("RMTProject/R/myRfunctions.R")

mu=1
sigma=0.5
d=5
S=50
N=100

dlist=seq(4,10, by=0.25)
mulist=seq(-3,3, by=0.25)
siglist=seq(0,1.5,by=0.1)

params=expand.grid(mulist,siglist,dlist)
k=nrow(params)
r=numeric(k)
l=r
r2=r
l2=r
r3=r
l3=l

for (i in 1:k)
{

    A=makeASymmetric(mu=params[i,1],sigma=params[i,2],d=params[i,3],S=S)#createA
    
  if(errorPD(A)==1){
    print("A is not positive definite")
  }else{
    X=makeX(mu=params[i,1],S=S,N=N,A=A)#createX
    
    r[i]=rhsD(S=S,N=N,X=X)#rhs for D in plogP wrt D
    l[i]=pDetD(sigma=params[i,2],d=params[i,3])#plogDet wrt D
    
    r2[i]=rhsmu(S=S,N=N,X=X)# rhs for mu in plogP wrt mu
    l2[i]=pDm(mu=params[i,1],d=params[i,3])#plogDet wrt mu
    
    r3[i]=rhssigma(S=S,N=N,X=X)
    l3[i]=pDetS(sigma=params[i,2],d=params[i,3])
    
    print(i)
  }
}

write.csv(c(r,l), "RMTProject/Data/Optimization/rhs vs lhs in plogP wrt d.txt")
write.csv(c(r2,l2), "RMTProject/Data/Optimization/rhs vs lhs in plogP wrt mu.txt")
write.csv(c(r3,l3), "RMTProject/Data/Optimization/rhs vs lhs in plogP wrt sigma.txt")


#############################
#Plots
##########

#d
plot(r,l, xlab="rhs in optimization for d", ylab="lhs in optimization for d", col="blue",
     title(main="Rhs vs Lhs in optimization of P wrt d"))
#mu
plot(r2,l2, xlab="rhs in optimization for mu", ylab="lhs in optimization for mu", col="blue",
     title(main="Rhs vs Lhs in optimization of P wrt d"))
#sigma
plot(r3,l3, xlab="rhs in optimization for sigma", ylab="lhs in optimization for sigma", 
     col="blue",title(main="Rhs vs Lhs in optimization of P wrt d"))


A=makeASymmetric(mu=mu,sigma=sigma,d=d,S=S)#createA

X=makeX(mu=mu,S=S,N=N,A=A)#createX

######################################
#for mu,sigma,d as defined at the top
#######################################

#mu
rhsmu(S=S,N=N,X=X)#rhs for mu in plogP wrt mu
pDm(mu=mu,d=d)#plogDet wrt mu

#d
rhsD(S=S,N=N,X=X)#rhs for D in plogP wrt D
pDetD(sigma=sigma,d=d)#plogDet wrt D

#sigma
rhssigma(S=S,N=N,X=X)#rhs for sigma in pLogP wrt sigma
pDetS(sigma=sigma,d=d)#plogDet wrt sigma
