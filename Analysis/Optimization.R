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


# for (i in 1:k)
# {
#   
#   A=makeASymmetric(mu=params[i,1],sigma=params[i,2],d=params[i,3],S=S)#createA
#   
#   if(errorPD(A)==1){
#     print("A is not positive definite")
#   }else{
#     X=try(makeX(mu=params[i,1],S=S,N=N,A=A))#createX
#     if(class(X)=='try-error'){
#       r[i]=l[i]=r2[i]=l2[i]=r3[i]=l3[i]=NA
#     }else{
#       
#       r[i]=rhsD(S=S,N=N,X=X)#rhs for D in plogP wrt D
#       l[i]=pDetD(sigma=params[i,2],d=params[i,3])#plogDet wrt D
#       
#       r2[i]=rhsmu(S=S,N=N,X=X)# rhs for mu in plogP wrt mu
#       l2[i]=pDm(mu=params[i,1],d=params[i,3])#plogDet wrt mu
#       
#       r3[i]=rhssigma(S=S,N=N,X=X)
#       l3[i]=pDetS(sigma=params[i,2],d=params[i,3])
#     }
#     print(i)
#   }
# }
# 
# write.csv(cbind(r1,l1), "RMTProject/Data/Optimization/rhs vs lhs in plogP wrt d2.txt")
# write.csv(cbind(r21,l21), "RMTProject/Data/Optimization/rhs vs lhs in plogP wrt mu2.txt")
# write.csv(cbind(r31,l31), "RMTProject/Data/Optimization/rhs vs lhs in plogP wrt sigma2.txt")
# 

#############################
#Plots
##########

# #d
# plot(r,l, xlab="rhs in optimization for d", ylab="lhs in optimization for d", col="blue",
#      title(main="Rhs vs Lhs in optimization of P wrt d"))
# #mu
# plot(r2,l2, xlab="rhs in optimization for mu", ylab="lhs in optimization for mu", col="blue",
#      title(main="Rhs vs Lhs in optimization of P wrt d"))
# #sigma
# plot(r3,l3, xlab="rhs in optimization for sigma", ylab="lhs in optimization for sigma", 
#      col="blue",title(main="Rhs vs Lhs in optimization of P wrt d"))


####################################################################################

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

#############################################
#solve for mu,sigma,d
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

n <- 20

params <- expand.grid(seq(sigma - 5, sigma + 5, length.out = n), 
                      seq(d - 5, d + 5, length.out = n), 
                      seq(mu - 5, mu + 5, length.out = n))
params <- as.matrix(params)

foo <- sapply(1:nrow(params), function(i) {
    Opfunct(c(params[i, ], mu), S, N, X)
})


plot(params[, 2], foo[3, ])
abline(h = 0)
abline(v = sol2$x[2])



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
# ep
# mu=params[376,1]
# d=params[376,3]
# sigma=params[376,2]
# params[376:377,]


#multiroot(f=model,S=S,N=N,X=X,start=c(0,1))

