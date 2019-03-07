mu=2
sigma=1.5
d=15
S=30
N=100
B=100

X=makeX(mu,sigma,d,S,N)
optim(par=c(mu,sigma,d),fn=.mc2opt,X=X,B=B, method = 'L-BFGS-B', 
      lower = c(-10, 0, 0), upper = c(10, 10, 100),
      control = list(fnscale = -1))
