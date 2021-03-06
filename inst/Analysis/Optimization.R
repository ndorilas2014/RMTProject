source("RMTProject/R/myRfunctions.R")

mu=1
sigma=0.5
d=5
S=60
N=100

dlist=seq(4,10, by=0.25)
mulist=seq(-3,3, by=0.25)
siglist=seq(0,1.5,by=0.1)

params=expand.grid(mulist,siglist,dlist)
k=nrow(params)
r=numeric(k)
l=r=l2=r2=l3=r3=l1=r1=l21=r21=l31=r31=r


for (i in 1:k)
{

  A=makeASymmetric(mu=params[i,1],sigma=params[i,2],d=params[i,3],S=S)#createA

  if(errorPD(A)==1){
    print("A is not positive definite")
  }else{
    X=try(makeX(mu=params[i,1],S=S,N=N,A=A))#createX
    if(class(X)=='try-error'){
      r[i]=l[i]=r2[i]=l2[i]=r3[i]=l3[i]=NA
    }else{

      r[i]=rhsD(S=S,N=N,X=X)#rhs for D in plogP wrt D
      l[i]=pDetD(sigma=params[i,2],d=params[i,3])#plogDet wrt D

      r2[i]=rhsmu(S=S,N=N,X=X)# rhs for mu in plogP wrt mu
      l2[i]=pDm(mu=params[i,1],d=params[i,3])#plogDet wrt mu

      r3[i]=rhssigma(S=S,N=N,X=X)
      l3[i]=pDetS(sigma=params[i,2],d=params[i,3])
    }
    print(i)
  }
}

# write.csv(cbind(r1,l1), "RMTProject/Data/Optimization/rhs vs lhs in plogP wrt d2.txt")
# write.csv(cbind(r21,l21), "RMTProject/Data/Optimization/rhs vs lhs in plogP wrt mu2.txt")
# write.csv(cbind(r31,l31), "RMTProject/Data/Optimization/rhs vs lhs in plogP wrt sigma2.txt")
# 

# N1=read.csv("RMTProject/Data/Optimization/rhs vs lhs in plogP wrt d.txt")
# N2=read.csv("RMTProject/Data/Optimization/rhs vs lhs in plogP wrt mu.txt")
# N3=read.csv("RMTProject/Data/Optimization/rhs vs lhs in plogP wrt sigma.txt")

#############################
#Plots
##########

# #d
# plot(N1[,2],N1[,3], xlab="rhs in optimization for d", ylab="lhs in optimization for d", col="blue",
#      title(main="Rhs vs Lhs in optimization of P wrt d"))
# 
# #mu
# plot(N2[,2],N2[,3], xlab="rhs in optimization for mu", ylab="lhs in optimization for mu", col="blue",
#      title(main="Rhs vs Lhs in optimization of P wrt mu"),xlim=c(0,5))
# 
# #sigma
# plot(N3[,2],N3[,3], xlab="rhs in optimization for sigma", ylab="lhs in optimization for sigma",
#      col="blue",title(main="Rhs vs Lhs in optimization of P wrt sigma"),xlim=c(0,0.15))


####################################################################################

A=makeASymmetric(mu=mu,sigma=sigma,d=d,S=S)#createA

X=makeX(mu=mu,S=S,N=N,A=A)#createX

######################################
#for mu,sigma,d as defined at the top
#######################################

#mu
rhsmu(S=S,N=N,X=X)#rhs for mu in plogP wrt mu
pDm(mu=2,d=d)#plogDet wrt mu

#d
rhsD(S=S,N=N,X=X)#rhs for D in plogP wrt D
pDetD(sigma=sigma,d=d)#plogDet wrt D

#sigma
rhssigma(S=S,N=N,X=X)#rhs for sigma in pLogP wrt sigma
pDetS(sigma=sigma,d=d)#plogDet wrt sigma

#############################################
#solve for mu,sigma,d using multiroot
sol=multiroot(f=Opfunct, S=S, N=N, X=X ,start=c(log(sigma), log(d), mu))


# solving using a different function from a different package that seems to work better for us
sol2 = nleqslv::nleqslv(c(log(sigma), log(d), mu), Opfunct, S = S, N = N, X = X, 
                        control = list(maxit = 300, allowSingular = TRUE))

# estimated sigma
exp(sol2$x[1])

# estimated d
exp(sol2$x[2])

# estimated mu
sol2$x[3]

# what we gave it
sigma
d
mu

# verifiy optimization in `sol2`

n <- 25

params <- expand.grid(seq(log(sigma) - 5, log(sigma) + 5, length.out = n), 
                      seq(log(d) - 5, log(d) + 5, length.out = n), 
                      seq(mu - 5, mu + 5, length.out = n))
params <- as.matrix(params)

foo <- sapply(1:nrow(params), function(i) {
    Opfunct(params[i, ], S, N, X)
})

##hERERER!!!

plot(exp(params[, 1]), foo[1, ],col="blue", xlab="values of sigma", ylab="optimization ouput", title(main="Optimization function vs. Sigma"))
 abline(h=0)
###########
plot(params[, 3], foo[3, ])
abline(h = 0)
abline(v = sol2$x[2])

###

#optimized mu,sigma, d
sigma1=sol$root[1]
d1=sol$root[2]
mu1=sol$roots[3]

################################################
##plot rhs/lhs vs min eg value =min(d-2sigma,d-mu)
w=numeric(nrow(params))
w2=w

for (i in 1:nrow(params))
{
  mu=params[i,1]
  sigma=params[i,2]
  d=params[i,3]
  w2[i]=min(d-2*sigma,d-mu)
  print(i)
}

plot(w,r/l)
q=r/l

##go through params, if d-2*sigma or d-mu<0 then save those parameter values
ep=numeric(nrow(params))
for (i in 1:nrow(params))
{
  mu=params[i,1]
  sigma=params[i,2]
  d=params[i,3]  
  if(q[i]>1.3 || is.na(q[i]))
  {
    ep[i]=i
  }else{
    ep[i]=0
  }
  print(i)
}


#generating combinations of mu,sigma, d


sigma=seq(0.01,1.5,length.out=10)
d=seq(4,10,length.out=10)
mu=seq(-3,3,length.out=10)
params=expand.grid(log(sigma),log(d), mu)
params=as.matrix(params)



#helper function to generate the data and run optimizatuin given vector of parameters p
newfoo=function(p){
  A=makeASymmetric(p[3],exp(p[1]),exp(p[2]),S)
  X=makeX(p[3], S,N,A)
  sol2 = nleqslv::nleqslv(c((p[1]), (p[2]), p[3]), Opfunct, 
                              S = S, N = N, X = X, 
                              control = list(maxit = 300, allowSingular = TRUE))
  return(sol2$x)
}
#getting the predicted value of mu, sigma, d
predictedP <- sapply(1:nrow(params), function(i) {
  print(i)
  output=try(newfoo(params[i,]))

  if(class(output)=="try-error"){
    return(rep(NA,3))
  }else{
    return(output)
  }
 
  
})

#plotting actual vs predicted d
plot(exp(params[,2]),exp(predictedP[2,]), xlab="actual value of d", 
     ylab="predicted value of d", col="red", title(main="Predicted vs. Actual d"))


foo3=matrix(numeric(nrow(params)),nrow=nrow(params),ncol=ncol(params))
###Plots for actual vs expected d
###plots partial of likelihood wrt mu,sigma, d
for(i in 1:nrow(params)){
foo3[i]=Opfunct(c(params[i,1],params[i,2],params[i,3]),S,N,X)
}

