source("R/myRfunctions.R")

load("EmpiricalData/feces_F4.Rdata")
load("EmpiricalData/feces_M3.Rdata")



#Normalizing data

a<-t(t(feces_F4)/ colSums(feces_F4))
#a<-t(t(feces_M3)/ colSums(feces_M3))
X<-(a-rowMeans(a))/sqrt(rowMeans(a^2)-rowMeans(a)^2)

mu=mean(X)
sigma=sd(X)
d=mean(diag(X))
S=nrow(X)
N=ncol(X)

out<-optim(par=c(mu,sigma,d),fn=.mc2opt,X=X,B=B, method = 'L-BFGS-B', 
           lower = c(-10, 0, 0), upper = c(10, 10, 100),
           control = list(fnscale = -1))
out$par
c(mu,sigma,d)
