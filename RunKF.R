
Run_KF<-function(par,old_par,osp_state,ostate,otime,osw_ind,opon_ind, opon_time,olik, nstate,ntime,nsw_ind, npon_ind,npon_time,acceptm)
{
#(par,old_par,nstate,osp_state,ostate,ntime,nsw_ind,olik,acceptm)
  nalpha = par[1]
  nbeta  = par[2]
  nrho   = par[3]
  nsigma = par[4]
  ntheta_x = par[5]
  ntheta_y = par[6]
  nBsigma = par[7]
  
  nswlamda = par[-(1:7)] 
  #nswlamda = par[8:9]

  alen = length( ntime )  

threshd = 1e-5 ## 0 or 1e-6 cause numerical error
if(nalpha>threshd&nrho>threshd&nsigma>threshd&nBsigma>threshd)
 { 


#Pinf initial value of P  or \Delta in paper
     # rsde<-sde(nai,nalpha,nbeta,nrho,nsigma)
     # A=rsde$A
     # Pinf= rsde$Pinf

  rsde<-sdeBM(nai,nalpha,nbeta,nrho,nsigma)
  A=rsde$A
  Pinf= rsde$Pinf


#Pinf initial value of P
    P_t = Pinf
    
    #a_y=1 ; a_x=1
    #m_yt =  c(0, rep(ay[1,a_y], nai+1))    ### z matrix 
    #m_xt =  c(0, rep(ax[1,a_x], nai+1))
    a_y <- min(which(!is.na(animaly[1, ])))
    a_x <- min(which(!is.na(animalx[1, ])))   #only need one coordinate but included for completeness
    m_yt =  c(0, rep(animaly[1,a_y], nai+1))    ### z matrix 
    m_xt =  c(0, rep(animalx[1,a_x], nai+1))
    

    Y_t=c(animaly[1,])
    ZY_t=Y_t[!is.na(Y_t)]  
    
    X_t=c(animalx[1,])
    ZX_t=X_t[!is.na(X_t)]
    
    #ZY_t=c(ay[1,])
    #ZX_t=c(ax[1,])
    
    ZY_t = ZY_t[-1]
    ZX_t = ZX_t[-1]
    
    v_yt = ZY_t - (H_t[[1]][-1,])%*%m_yt
    v_xt = ZX_t - (H_t[[1]][-1,])%*%m_xt
    
    
    S_t = H_t[[1]][-1,]%*%P_t%*%t(H_t[[1]][-1,])
    detst = det(S_t)
    solst = solve(S_t)
    
    k_t = P_t%*%t(H_t[[1]][-1,])%*%solst
    
    m_yt = m_yt + k_t%*%(v_yt)
    m_xt = m_xt + k_t%*%(v_xt)
    P_t = P_t - k_t%*%S_t%*%t(k_t)
  

     sp_ind = 1
     # m_yt = matrix(c(ntheta_y,black[2,1],ay[sp_ind,]),ncol=1)
     # m_xt = matrix(c(ntheta_x,black[1,1],ax[sp_ind,]),ncol=1)
      lik = 0

    for(ia in 2:(alen) )
   {
      dt = ntime[ia] - ntime[ia-1]

      # rsol=solvesde(A,Pinf,dt)
      rsol = solvesdeBM(A,dt,nalpha,nbeta,nrho,nsigma)

      F = matrix(rsol$F,ncol= dim(rsol$F)[1] )
      Casigma = rsol$Casigma
      
      oF = F
      oCasigma = Casigma


      ouindex = which(c(nstate[ia-1,])==1)+2 ## ss[1,] is the starting state of movement at time 2
      bmindex = which(c(nstate[ia-1,])==2)+2

      if(length(ouindex)>0)
         {
            F[ouindex,] = oF[ouindex,]
            Casigma[ouindex,]=as.matrix(oCasigma[ouindex,]) 
            Casigma[,ouindex]=as.matrix(oCasigma[,ouindex])
         }
      
       if(length(bmindex)>0)
         {
            if(length(bmindex)>1)
            {
              ## set F row and collumn to zero for BM animal
              F[bmindex,] = 0
              diag(F[bmindex,bmindex])=1
              ## set Casigma row and collumn to zero for BM animal
              Casigma[bmindex,] = 0
              Casigma[,bmindex] = 0
              diag(Casigma[bmindex,bmindex]) = nBsigma^2
            } else if(length(bmindex)==1)
            {
              ## set F row and collumn to zero for BM animal
              F[bmindex,] = 0
              F[bmindex,bmindex] = 1
              ## set Casigma row and collumn to zero for BM animal
              Casigma[bmindex,] = 0
              Casigma[,bmindex] = 0
              Casigma[bmindex,bmindex] = nBsigma^2
            }
          }
     
    if( nsw_ind[ia] == 'SW')
    #if( any( is.na(animaly[ia,]) ) )
     {
   ## predict step
   ## predict new mu based on new A and new lambda mu_t = m_0 e^A + (I - e^A)*theta  
       m_yt = F%*%m_yt 
       m_xt = F%*%m_xt 
   ## predict new var based on new A and new lambda p_t = e^A p_t-1 e^A' + Q   nCsigma =Q
       P_t = F%*%P_t%*%t(F) + Casigma

     } 
     else if( nsw_ind[ia] == 'SP')
     {            
   ## predict step
   ## predict new mu based on new A and new lambda mu_t = m_0 e^A + (I - e^A)*theta  
       m_yt = F%*%m_yt 
       m_xt = F%*%m_xt 
   ## predict new var based on new A and new lambda p_t = e^A p_t-1 e^A' + Q   nCsigma =Q
       P_t = F%*%P_t%*%t(F) + Casigma
   ## updating step
       #Y_t=c(animaly[ia,]) 
       #X_t=c(animalx[ia,])
       sp_ind = sp_ind+1

       rY_t=c(ay[sp_ind,])
       Y_t=rY_t[!is.na(rY_t)] 

       rX_t=c(ax[sp_ind,])
       X_t=rX_t[!is.na(rX_t)]

       v_yt = Y_t - H_t[[sp_ind]]%*%m_yt
       v_xt = X_t - H_t[[sp_ind]]%*%m_xt
       S_t = H_t[[sp_ind]]%*%P_t%*%t(H_t[[sp_ind]])

       detst = det(S_t)
       solst = solve(S_t)

       logde = tryCatch({ log(detst)
        }, warning = function(war)
        { 
          print(paste("log det: ",war))
          return( logde = -Inf )
        },
         error = function(err) 
        {# error handler picks up where error was generated
        },finally = {  } )

        if(logde == -Inf)
        {
          lik = -Inf
          break
        }

        lik1 = as.matrix(- 1/2*nai*log(2*pi) - 1/2*logde -1/2*t(v_yt)%*%solst%*%v_yt) #p(y_t|y_1:t-1,theta)
        lik2= as.matrix(- 1/2*nai*log(2*pi) - 1/2*logde -1/2*t(v_xt)%*%solst%*%v_xt) #p(y_t|y_1:t-1,theta)
        lik = lik + lik1 + lik2                                           # p(y_1:t|theta)

        k_t = P_t%*%t(H_t[[sp_ind]])%*%solst 
        m_yt = m_yt + k_t%*%(v_yt)
        m_xt = m_xt + k_t%*%(v_xt)
        P_t = P_t - k_t%*%S_t%*%t(k_t)

     }

  }

## if we only update diffusion par, we leave switching probability. if we update state, we include switching probability
     #if( all( (par-old_par )!=0)  )
     #{
     # sw = osw
     # sigHR <- exp(lik- olik)
     #} else {
      # sw = 0
      # for(jjj in 1:nai)
      # {
      # sw=sw+sum(log(jp[ cbind(nstate[-alen,jjj],nstate[-1,jjj]) ]) ) 
      # }
      #  sigHR <- exp(sw+lik- olik-osw)
     #}
    sigHR <- exp(lik - olik )  
    
  } 
   else
  {
     sigHR<- 0
  }

  if(runif(1) < sigHR[1]) 
   {
       alpha = nalpha
       beta = nbeta
       rho = nrho
       sigma = nsigma
       theta = c(ntheta_x,ntheta_y)
       Bsigma = nBsigma
       swlamda = nswlamda
       
       acceptm = acceptm+1
       olik = lik 
       #olnlik_lam = nlnlik_lam
       #osw = sw
       sp_state = nstate[which(nsw_ind=='SP'),]
       ostate = nstate
       otime = ntime
       osw_ind = nsw_ind
       opon_ind  = npon_ind 
       opon_time  = npon_time 
   }
   else 
   {
       alpha = old_par[1]
       beta = old_par[2]
       rho = old_par[3]
       sigma = old_par[4]
       theta = c(old_par[5],old_par[6])
       Bsigma = old_par[7]
       swlamda = old_par[-(1:7)] 
       #swlamda = old_par[8:9]

       sp_state = osp_state
     }

  return( list("par"=c(alpha,beta,rho,sigma,theta,Bsigma,swlamda),"olik"=olik,"accept"=acceptm,"SPstate"=sp_state,"ostate"=ostate,"otime"=otime,"osw_ind"=osw_ind,"opon_ind"=opon_ind,"opon_time"=opon_time) )

}

