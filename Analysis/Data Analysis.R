##Analysis



#Data Analysis

#Data for A and X
source("Data/A&X_data.R")

#Functions
source("R/myRfunctions.R")

A=makeA_BV(mu=mu, sigma=sigma, rhoB=rhoB, d=d, S=S)

Eg=eigen(A)
minEg=min(Re(Eg$values))
maxEg=max(Re(Eg$values))
#If the real part of the right most eigenvalue is nonnegative then the matrix is unstable
if(maxEg>=0)
{
  print("Error: The matrix is not diagonally stable")
} else
  
{
  ###Variables
  ###
  T=500
  dt=0.01 #time steps
  y=matrix(0,nrow=S, ncol=T) #population
  y[,1]=30 #initial condition
  
  
  for (j in 2:T)
  {
    for (i in 1:S)
    {
      noise=rnorm(1)*10 #adding random noise variable
      y<-TS(A, dt, y, xi=noise, i, j)
    }
  }
  
  
  #matplot((c(1:T)*dt),t(y), lty=1, col=rainbow(S), type="l") #time series plot
  
  
  ################################
  ##STATISTICAL PROPERTIES OF Y
  ################################
  
  mbs=rnorm(T)
  ms=rnorm(S)
  vbs=rnorm(T)
  vs=rnorm(S)
  cvs=matrix(rnorm(S*S), S)
  
  
  ##mean & variation between species at time i
  ##
  for(i in (1:T))
  {
    mbs[i]=mean(y[,i]) #mean
    vbs[i]=var(y[,i]) #variation
  }
  #plot((c(1:T))*dt, mbs, type="l", col="blue")#plot the mean population between species at each time step
  #plot((c(1:T))*dt, vbs, type="l", col="blue")#plot the variance of pop size between species at each time step
  
  
  ##mean & variation of species i across time
  ##
  for (i in (1:S))
  {
    ms[i]=mean(y[i,])
    vs[i]=var(y[i,])
  }
  #plot(c(1:S), ms, type="l", col="red") #plot the mean population of each ind. species across time
  #plot(c(1:S), vs, type="l", col="red") #plot the variation in pop. of each ind. species across time
  
  
  
  for (i in (1:S))
  {
    for (j in (1:S))
    {
      cvs[i,j]=cor(y[i,], y[j,])#correlation between species i and the rest of the species(1-S)
    }
    # for(k in (1:T))
    # {
    #   
    # }
  }
  
  
  matplot(c(1:S), t(cvs), type="l", col=rainbow(S), lty=1) # plot corr between species i and rest
  #plot(c(1:S), cvs[1,], type = "l", col="red")
  #lines(c(1:S), cvs[2, ], type ="l", col="blue")
  
}