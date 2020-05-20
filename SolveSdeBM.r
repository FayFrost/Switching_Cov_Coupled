
solvesdeBM<-function(A,dt,alpha,beta,rho,sigma)
{
  # F = as.matrix(expm::expm(A*dt) )
  
  
  fe <- exp((-alpha)*dt)
  
  F = diag(nai+2)
  F = F*fe 
  F[3:(nai+2), 2] = 1 - fe
  F[1,1]=1
  F[2,2]=1
        
###  Var  = delta - e^At*delta*e^A't     Q variance   or  \Xi in paper
     # csigma=Pinf-F%*%Pinf%*%t(F)
     # Casigma=as.matrix(csigma)   

 Aid= diag(nai+2)

###### create  the stationary matrix   lambda matrix 
  ELF =(rho^2)*dt - ((rho^2)/alpha)*(1-exp((-alpha)*dt))
  EF  =((sigma^2/(2*alpha))*(1-exp(-2*alpha*dt))) + (((rho^2)/(2*alpha))*(2*alpha*dt - 3)) + (2*(exp(-alpha*dt))*rho^2)/alpha - (exp(-2*alpha*dt)*rho^2)/(2*alpha)
  EL = (rho^2)*dt 

  EFF = EF - (((sigma^2/(2*alpha))*(1-exp(-2*alpha*dt)))) 
  

  lamda =Aid*EF
  
  lamda[2:(nai+2),2:(nai+2)] = lamda[2:(nai+2),2:(nai+2)]+EFF - EFF*Aid[2:(nai+2),2:(nai+2)]
  lamda[2, 2:(nai+2)] = ELF
  lamda[2:(nai+2), 2] = ELF
  
  lamda[2,2] = EL


return(list(F=F,Casigma=lamda) )
}
