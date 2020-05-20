
library(mvtnorm)
library(Matrix)
set.seed(3)



alpha <- 0.5 #0.8
beta  <- 1e-6 
rho   <- 3  #5.2
sigma <- 1.0   # 2.7
Bsigma <- 0.5

ho = 0.5   ##rate parameter of logistic BM to OU
hb = 0.5  ##rate parameter of logistic OU to BM

tou = 20  ##most likely time of day to switch to OU
tb =   8 ##most likely time of day to switch to  BM 



covariates = TRUE
if(covariates){
  looklamda <- function(t)
  {
    #cat("looklambda called with",t,"\n")
    lam12 <-(ho/2)*(1+cos((t-tou)*(2*pi/24)))
    lam21 <-(hb/2)*(1+cos((t-tb)*(2*pi/24)))
    return(c(lam12, lam21))
  }
}else{
  looklamda <- function(t){
    swlamda = c(0.1, 0.4)     ## fixed switching rates 
    return(swlamda)
  }
}


#delata T, the sampling time or the grid time inverval  
# Nsamp how many samples we want
deltat=2  
dt = deltat ## used for inference code in other files
Nsamp = 50
samp = 1   

# how many animals
n = 6


Kappa = n*max(ho, hb)

#initial black dot and animal location
init_x = c( 1, 1.5, 1.5, -1, 1.7, 1.7, -1.2)
init_y = c( 1, 1.5, 1,   -1, 1.7, 1.2, -1.2)

# initial states
state<-rep(c(1),n+1)    # state 1 attract to black dot         state 2 unattract to black dot (brownian motion)
# first element of state is useless
Px = init_x
Py = init_y

swt_ls <- list( matrix(c(0)),matrix(c(0)),matrix(c(0)),matrix(c(0)),matrix(c(0)),matrix(c(0)) ) # was 5 matrix elements for 5 animals

sws_ls <- swt_ls 
sw_all <- list(matrix(c(0)), matrix(c(0)),matrix(c(0)))

sampl <- list( matrix(c(0)) )

sw_ind <- c('SP')
all_t <- c(0)
all_s <- state
all_x <- c(Px)
all_y <- c(Py)
swt_old =0

# A matrix   coefficient matrix of SDE
Aid= diag(n+1)
A=Aid*(-alpha)
A[1:n+1,1] = alpha
A[1,1]= -beta

covA=rho^2*alpha/(2*beta*(alpha+beta))
varA=sigma^2/(2*alpha)+rho^2*alpha/(2*beta*(alpha+beta))
lamda =Aid*varA  
lamda[1:(n+1),1:(n+1)] =lamda[1:(n+1),1:(n+1)]+covA-covA*Aid[1:(n+1),1:(n+1)]
lamda[1,1] = rho^2/(beta*2)

## record states at sampling point for all animal
state_sp = array(c(1),dim=c(Nsamp,n))  
location_sp = array( (0),dim = c(Nsamp,n+1,2) )


#rest T ,how far from current time to the grid time point
resT= deltat

ouindex=which(state==1); n_ou <- sum(state[-1]==1)    ##count number of FOLLOWERS in OU state, to be used in P_ac_sw
bmindex=which(state==2); n_bm = length(bmindex)

pon_probAll = c()
## propose first switch time
swt = 0
npsw = 0
repeat{
  po_t = rexp(1, Kappa)
  swt = swt + po_t            #potential switching time 
  swlamda_now <- looklamda(swt) ## to use correct switching rates for the time point
  #cat("swlambda_now", swlamda_now, "time", pswt + up_start,  "\n")
  P_ac_sw =( n_ou*swlamda_now[1] + n_bm*swlamda_now[2] )/Kappa # prob of actual switch
  
  ac_sw = sample( c(0,1), 1, FALSE, prob=c(1-P_ac_sw,P_ac_sw) ) # test for actual switch
  npsw=npsw+1
  if(ac_sw==1) break else
    pon_probAll <- c(pon_probAll,1-P_ac_sw/Kappa) 
}


# swt=rexp(1,n*looklamda(0)[1])  #switch time

resT = swt - deltat
move2 = swt
move1 = deltat

qq = 0
tt = 0


tts = swt   ###first switch 
tto = deltat

cumT <- 0

repeat
{
  ##case 1, switch animal state
  
  if(resT<0)     ## switching happen first , before sampling
  {
    
    
    cumT = cumT + tts
    tto = tto - tts
    #cat("switching happen first############# \n")
    ### move to switch point
    t = move2  
    #cat("movement at time",t,"\n")
    ### make the move
    ouindex=which(state==1); n_ou <- sum(state[-1]==1)    ##count number of FOLLOWERS in OU state, to be used in P_ac_sw
    bmindex=which(state==2); n_bm = length(bmindex)
    
    if(length(ouindex)>1)  
    {      ## move black dot and animal with OU
      csigma=lamda[ouindex,ouindex]-expm(A[ouindex,ouindex]*t)%*%lamda[ouindex,ouindex]%*%expm(t(A[ouindex,ouindex])*t)   # var  c- ACA
      csigma=trunc(csigma*10^5)/10^5
      Casigma=as.matrix(csigma)     
      
      mux<-expm(A[ouindex,ouindex]*t)%*%Px[ouindex]                           # mean
      muy<-expm(A[ouindex,ouindex]*t)%*%Py[ouindex]
      cat("animal",ouindex," move OU \n")
      
      Px[ouindex]=rmvnorm(1,as.matrix(mux),Casigma,"chol")
      Py[ouindex]=rmvnorm(1,as.matrix(muy),Casigma,"chol")   
    }
    
    if(length(bmindex)>0)
    {     ## move animal with BM
      #cat("animal",bmindex,"move BM \n")		   
      bmcov=diag(Bsigma^2,length(bmindex))*t
      
      Px[bmindex]=rmvnorm(1,Px[bmindex],bmcov,"chol")
      Py[bmindex]=rmvnorm(1,Py[bmindex],bmcov,"chol")
    }
    
    
    ## pick the switch animal 
    curr_lamda <- looklamda(cumT)
    de=sum(curr_lamda[state])
    # SwitchAnimal =sample(c(2:(n+1) ),1,FALSE,prob=c( swlamda[state[2]]/de,swlamda[state[3]]/de,swlamda[state[4]]/de,swlamda[state[5]]/de,swlamda[state[6]]/de,swlamda[state[7]]/de ))
    SwitchAnimal =sample( c(2:(n+1) ),1,FALSE,prob=c(curr_lamda[ state[ 2:(n+1) ] ]/de) )
    
    
    ## switch animal state 
    state[SwitchAnimal] = 3 - state[SwitchAnimal]
    
    ## demo state and switch time
    swt_ls[[SwitchAnimal-1]] <- cbind( swt_ls[[SwitchAnimal-1]], swt ) 
    sws_ls[[SwitchAnimal-1]] <- cbind( sws_ls[[SwitchAnimal-1]], state[SwitchAnimal] ) 
    sw_all[[1]]<- cbind( sw_all[[1]],swt )
    sw_all[[2]]<- cbind( sw_all[[2]],state[SwitchAnimal] )
    sw_all[[3]]<- cbind( sw_all[[3]], (SwitchAnimal-1) )
    cat(tt,"iteration ", "switch", swt ,"\n")
    
    sw_ind= rbind(sw_ind,'SW')
    all_s = rbind(all_s,state)
    all_t = cbind(all_t,swt+ sum(swt_old) )
    all_x = rbind( all_x, rep(NA,n+1) )
    all_y = rbind( all_y, rep(NA,n+1) )
    swt_old = c(swt_old,swt)
    
    ## update index   
    
    ouindex=which(state==1); n_ou <- sum(state[-1]==1)    ##count number of FOLLOWERS in OU state, to be used in P_ac_sw
    bmindex=which(state==2); n_bm = length(bmindex)
    ## propse next switch 
    # cat("ouindex",ouindex,"bmindex",bmindex,"sum","\n")
    
   # swt=rexp(1, sum(length(ouindex)*looklamda(cumT)[1],length(bmindex)*looklamda(cumT)[2]) ) 
   
    swt = 0      
    npsw = 0
    repeat{
      po_t = rexp(1, Kappa)
      swt = swt + po_t
      swlamda_now <- looklamda(cumT + swt) 
      
      P_ac_sw =( n_ou*swlamda_now[1] + n_bm*swlamda_now[2] )/Kappa
      
      ac_sw = sample( c(0,1), 1, FALSE, prob=c(1-P_ac_sw,P_ac_sw) )
      npsw=npsw+1
      
      cat("swlamda_now", swlamda_now, "time", swt + cumT, "npsw", npsw,  "\n")
      if(ac_sw==1) break else
        pon_probAll <- c(pon_probAll,1-P_ac_sw/Kappa)
    }
    
    
  
    
    move1 = abs(resT)
    move2 = swt
    
    
    tts = swt
    resT = swt - abs(resT)
    
    
    #cat(resT,"rest time to next sampling point \n")
    #cat("state",state[2]," ",state[3]," ",state[4]," ",state[5]," ",state[6]," ",state[7]," \n \n")
    # cat(tt,"iteration ", "switch", swt , "cumT", cumT, "\n") 
    
    cat( "cumT", cumT, "\n") 
  }
  
  ###case 2, do not switch animal state
  if(resT > 0)            #sampling happan first  resT>0
  {
    
    
    # cat("sampling happen first############# \n")
    ### move to sampling point
    t = move1  
    # cat("movement at time",t,"\n")
    ### make move
    ouindex=which(state==1); n_ou <- sum(state[-1]==1)    ##count number of FOLLOWERS in OU state, to be used in P_ac_sw
    bmindex=which(state==2); n_bm = length(bmindex)
    
    if(length(ouindex)>1)  
    {
      csigma=lamda[ouindex,ouindex]-expm(A[ouindex,ouindex]*t)%*%lamda[ouindex,ouindex]%*%expm(t(A[ouindex,ouindex])*t)   # var  c- ACA
      csigma=trunc(csigma*10^5)/10^5
      Casigma=as.matrix(csigma)     
      
      mux<-expm(A[ouindex,ouindex]*t)%*%Px[ouindex]                           # mean
      muy<-expm(A[ouindex,ouindex]*t)%*%Py[ouindex]
      #cat("animal",ouindex,"move OU \n")
      
      Px[ouindex]=rmvnorm(1,as.matrix(mux),Casigma,"chol")
      Py[ouindex]=rmvnorm(1,as.matrix(muy),Casigma,"chol")
    }
    
    if(length(bmindex)>0)
    {
      #cat("animal",bmindex,"move BM \n")
      bmcov=diag(Bsigma^2,length(bmindex))*t
      
      Px[bmindex]=rmvnorm(1,Px[bmindex],bmcov,"chol")
      Py[bmindex]=rmvnorm(1,Py[bmindex],bmcov,"chol")
    }
    
    ### recording location and state information# 
    if(samp == Nsamp)
    {
      break
    } else
    {
      # save location and state at observation points
      location_sp[samp+1,,] = cbind(Px,Py)
      state_sp[(samp+1),(1:n)] = state[2:(n+1)]
      
      # demo state and sampling time
      demot = samp*2
     # cat(tt,"iteration sample at", demot ,  "\n")  
      
      samp = samp+1
      sampl = cbind( sampl, demot)
      
      sw_ind= rbind(sw_ind,'SP')
      all_s=  rbind(all_s,state)
      all_t = cbind(all_t,demot)
      all_x = rbind( all_x, Px )
      all_y = rbind( all_y, Py )
      #cat("state",state[2]," ",state[3]," ",state[4]," ",state[5]," ",state[6]," ",state[7]," \n \n")
    }
    
    move2 = abs(resT)
    move1 = deltat
    
    resT = resT - deltat
    cumT = cumT + tto
    tts = tts - tto
    tto = deltat 
    
    
    cat( "cumT", cumT, "\n") 
    #cat(resT,"rest time to next sampling point \n")
  }
  
  tt=tt+1
}




ss<- state_sp #t( rbind(s1,s2,s3,s4,s5,s6) )

ax <- location_sp[,-1,1] #cbind(a1[1,], a2[1,], a3[1,], a4[1,], a5[1,], a6[1,])
ay <- location_sp[,-1,2] #cbind(a1[2,], a2[2,], a3[2,], a4[2,], a5[2,], a6[2,])  

animalx = ax
animaly = ay

theta_x = 0
theta_y = 0
theta=c(theta_x,theta_y)


animation_plot<-function( break_seconds)
{
  n = dim(location_sp)[2]-1
  colour_l= c('black','red','blue','green','red','blue','green', 'orange')
  pl_state = state_sp
  pl_state[which(pl_state>1.5)] = 17      # BM
  pl_state[which(pl_state<1.5)] = 19      # OU
  
  plot(0,0,xlim=c( min(location_sp[,,1])-2, max(location_sp[,,2])+2 ),ylim=c( min(location_sp[,,2])-2, max(location_sp[,,2])+2 ),xlab="x-dimension",ylab="y-dimension")
  
  for(qq in 1:Nsamp)
  {
    
    for(ia in 1:n)
    {
      lines(location_sp[qq,ia,1],location_sp[qq,ia,2],type="p",pch=pl_state[qq,ia],col=colour_l[ia])
    }
    
    Sys.sleep(break_seconds)
    
    for(ia in 1:n)
    {
      lines(location_sp[qq,ia,1],location_sp[qq,ia,2],type="p",pch=pl_state[qq,ia],col="white",cex=2)
    }
    
  }
}

