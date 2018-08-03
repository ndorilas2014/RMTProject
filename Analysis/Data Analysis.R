##Analysis



#Data Analysis

#Data for A and X
source("RMTProject/Data/A&X_data.R")

#Functions
source("RMTProject/R/myRfunctions.R")
S=50
A=makeA_BV(mu=mu, sigma=sigma, rhoB=rhoB, d=d, S=S)
A2=makeASymmetric(mu=mu,sigma=sigma,d=-d, S=S)

#If A not negative definite then the matrix is unstable
if(errorND(A)==1)
{
  print("Error: The matrix is not diagonally stable")
} else
  
{
  ###Variables
  ###
  T=500
  dt=0.01 #time steps
  y=matrix(0,nrow=S, ncol=T) #population
  y[,1]=15 #initial condition
  
  
  for (j in 2:T)
  {
    for (i in 1:S)
    {
      noise=rnorm(1)*10 #adding random noise variable
      y<-TS(A, dt, y, xi=noise, i, j)
    }
  }
  
  
  matplot((c(1:T)),t(y),xlab="time(t)", ylab="population of species at time t", 
          main="Population of S Species vs. Time",lty=1, col=rainbow(S), type="l") #time series plot
  
  
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
  #plot((c(1:T)), mbs, type="l", xlab="time(t)",ylab="mean population", 
  #main="The Mean Population of Individuals at Time t", col="blue")#plot the mean population between species at each time step
  # plot((c(1:T)), vbs, type="l", xlab="time(t)",ylab="variance in population size", 
  #      main="Variance in Population Size Across Species at Time t", col="blue")#plot the variance of pop size between species at each time step
  # 
  
  ##mean & variation of species i across time
  ##
  for (i in (1:S))
  {
    ms[i]=mean(y[i,])
    vs[i]=var(y[i,])
  }
  plot(c(1:S), ms,  xlab="species #",ylab="mean across time",
  main="Mean in Population Size of Each Species Across Time", type="l", col="red") #plot the mean population of each ind. species across time
  plot(c(1:S), vs,  xlab="species #",ylab="variation across time",
  main="Variation in Population Size of Each Species Across Time", type="l", col="red") #plot the variation in pop. of each ind. species across time
  

  
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
  
  
  #matplot(c(1:S), t(cvs), type="l", col=rainbow(S), lty=1) # plot corr between species i and rest
  plot(c(1:S), cvs[1,], type = "l", col="red")
  lines(c(1:S), cvs[2, ], type ="l", col="blue")
  
}
