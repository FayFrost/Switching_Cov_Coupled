

sdeBM<-function(nai,alpha,beta,rho,sigma)
{
###  SDE   
##covariance matrix d x = A dx +  matrix(rho, sigma) dB    Casigma is cov matrix
#A= matrix(c(-beta,alpha,alpha,alpha,0,-alpha,0,0,0,0,-alpha,0,0,0,0,-alpha),nrow=4,ncol=4)
   
    #lamda = matrix(c(0),nrow=nai+2,ncol=nai+2)

#### initial matrix A     or A in paper
    Aid= diag(nai+2)
    Aid[1,1] = 0
    A=Aid*(-alpha)
    A[2:(nai+2),2] = alpha
    A[2,2]= -beta 
    A[2,1]= beta

#### create delta matrix    the stationary matrix   lambda matrix   or \Delta in paper
    #covA=rho^2*alpha/(2*beta*(alpha+beta))
    #varA=sigma^2/(2*alpha)+rho^2*alpha/(2*beta*(alpha+beta))
    #lamda =Aid*varA  
    #lamda[2:(nai+2),2:(nai+2)] =lamda[2:(nai+2),2:(nai+2)]+covA-covA*Aid[2:(nai+2),2:(nai+2)]
    #lamda[2,2] = rho^2/(beta*2)


######## Create Initial P 
    a_y <- min(which(!is.na(ay[1, ])))
    a_x <- min(which(!is.na(ax[1, ])))   #only need one coordinate but included for completeness
    

    P = diag(nai+2)
    P[2:(nai+2),2:(nai+2)] = (sigma^2)/(2*alpha) 
    P = P + ((sigma^2)/(alpha)*diag(nai+2)) - ((sigma^2)/(2*alpha)*diag(nai+2)) 
    P[1,] = 0
    P[a_y+2, ] = P[, a_y+2] = 0 
    P[2, 2] = (sigma^2 + 3*rho^2)/(2*alpha)


 return(list(A=A,Pinf=P))
}



