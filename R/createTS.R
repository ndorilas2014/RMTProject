#' @title Species Abundance Time Series
#' 
#' @description createTS creates time series data either through independent draws from a 
#' multivariate normal distribution with a symmetric interaction matrix A, using a model 
#' for species abundance(see details) with a symmetric interaction matrix A and using 
#' the same model but with a nonsymmetric A 
#' 
#' @param mu the mean of the elements in the interaction matrix A
#' @param sigma the deviatin of the elements in the interaction matrix A
#' @param d the diagonal entries of the elemnts in the interaction matrix A
#' @param rho the correlation (assuming A is created from a bivariate distribution)
#' @param dt the time steps that the species abundance model will be solved for 
#' @param y0 vector holding the initial condition for the species abundance(time=0). 
#' @param ndata the number of data points you would like to have 
#' @param t_err the deviation in the noise in the time series model(assuming noise~N(0,t_err))
#' @param opt has several options, see details below 
#' 
#' @details The opt function has three possible values: 
#' 1.'symmetricA_mvnorm' which allows for the creation of time series data for a symmetric 
#' interaction matrix and species data created from a multivariate normal distribution. 
#' 
#' 2. 'symmetricA_TS' which allows for the creation of time series data with a 
#' symmetric interaction  matrix and species abundance data from the model 
#' (1) y_i(t+dt)=y_i(t)+dt(sum_i,j(sum_t(x_i A x_J))). 
#' 
#' 3. 'bivariateA_TS' which uses a non symmetric interactioni matrix drawn from a bivariate 
#' distribution, and the time series model (1) to create the time series data.
#'  
#' @return A matrix where rows are species and columns are the time points
#' @export 



createTS<-function(mu,sigma,d,sizeA,rho=NULL,dt=NULL,y0=NULL,ndata,t_err=NULL,opt=c('symmetricA_mvnormal','symmetricA_TS', 'bivariateA_TS')){
    opt=match.arg(opt,c('symmetricA_mvnormal','symmetricA_TS', 'bivariateA_TS'))
    if(opt=='symmetricA_mvnormal')
    {
        return(.xtseriesmvnormal(mu,sigma,d,S=sizeA,N=ndata))
    }else if(opt=='symmetricA_TS'){
        A=makeASymmetric(mu,sigma,d,S)
        return(.xtseries(A,dt,y0,N=ndata,t_err))
    }else{
        A=makeA_BV(mu,sigma,rho,d,S)
        return(.xtseries(A,dt,y0,N=ndata,t_err))
    }
}





.xtseriesmvnormal<-function(mu, sigma, d, S, N)
{
    # browser()
    mulist=rep(mu, S)#Generating means for X
    A=makeASymmetric(mu,sigma,d,S)
    
    X=(t(rmvnorm(N, (mulist/S), inv(A), method="svd"))) #Generates X/data from a mv distribution
    #N=ncol(X) #for when X is given not created from mv distribution
    
    return(X)
}


.xtseries<-function(A,dt,y0,N,sigma)
{
    S=nrow(A)
    
    y=matrix(0,nrow=S, ncol=N) #population
    y[,1]=y0 #initial condition
    x1=rnorm(n=S, mean=0, sd=sigma)
    x2=matrix(rnorm(n=N*S, mean=0, sd=sigma),S,N)
    
    
    for(i in 1:(N-1))
    {
        
        y[,i+1]=y[,i]+dt*(A%*%y[,i]+x2[,i+1])
        
    }
    return(y)
    
}