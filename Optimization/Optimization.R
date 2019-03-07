mu=1
sigma=0.5
d=10
S=10
N=10
B=100

X=makeX(mu,sigma,d,S,N)
optim(par=c(mu,sigma,d),fn=.mc2opt,X=X,B=B, method = 'L-BFGS-B', 
      lower = c(-10, 0, 0), upper = c(10, 10, 100),
      control = list(fnscale = -1))
