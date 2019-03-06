mu=1
sigma=0.5
d=10
S=10
N=10
B=100

X=makeX(mu,sigma,d,S,N)
optim(par=c(mu,sigma,d),fn=MonteCarlo,X=X,B=B)
