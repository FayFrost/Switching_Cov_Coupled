
#install.packages("lubridate")
library(lubridate)

library(mvtnorm)
library(Matrix)
set.seed(3809)

source("SdeBM.r")
source("SolveSdeBM.r")
source('RunKF.r')


source('data_sim_cov.r')
#source('Data_real.r')

mxsamp = max(all_t)
alen = max(dim(ax))
all_state <- all_s
all_time <- as.matrix( all_t )
nai = 6 # 5 animals


covariates = FALSE # using covariates?



# initialise paramters 
alpha <- 0.5 #0.8
beta  <- 1e-6 
rho   <- 3 #5.2
sigma <- 1  # 2.7
Bsigma <- 0.5

## proposal SDs
sf = 2
prop.alpha = 0.01 * sf
prop.beta  = 0.0 #0.5 #0.1  * sf
prop.rho   = 0.1 * sf
prop.sigma = 0.1  * sf
prop.theta = 0.1  * sf
prop.Bsigma = 0.1 * sf

swlamda <- c(0.1, 0.4)
#initialise switchign rate paramters
ho <- 0.5 
hb <- 0.5

tou <- 20
tb <- 8

#prop sd
prop.ho <-0.005
prop.hb <- 0.005
  
prop.tou <- 0.0
prop.tb <- 0.0

#initialise switching rates if fixed, prop sd are in lookup lamda function
lamdaOB <- 0.1
lamdaBO <- 0.4



if(covariates){
  swlamda <- c(ho, hb, tou, tb)
  sw_prop_SD <- c(prop.ho, prop.hb, prop.tou, prop.tb)
  lookup_lamda <- function(swlamda, t)
  {
    lam12 <- (ho/2)*(1+cos((t-tou)*(2*pi/24)))
    lam21 <- (hb/2)*(1+cos((t-tb)*(2*pi/24)))
    return(cbind(lam12, lam21))
  }
}else{
  swlamda <- c(lamdaOB ,lamdaBO)
  sw_prop_SD <- c(0.05, 0.1)
  lookup_lamda <- function(swlamda, t){
       return(cbind(swlamda[1], swlamda[2]))
  }
}
## we fix Kappa, and also if the swiching rate from OU to BM is too small (close to zero), then it would be 
## hard to find a actual switching point. Because for many cases, animal would do OU and the actual switchig 
## probability become nai*lambda_ou_to_bm/kappa, which is very small. Therefore we need set low bound and high
## bound for lambdas


Kappa =  3.5
sw = 0
osw = -10^4


################## initialise accept number ################################################## 
lamaccept = 0
alaccept  = 0
#beaccept  = 0
sigaccept = 0
rhaccept  = 0
#theaccept = 0
Bsigaccept = 0
staccept = 0

dat <- strftime(Sys.time(),format="%Y%b%d_%H-%M")
#filetheta3 <- paste("blacBGtheta",dat,".txt", sep = "")
filetheta4 <- paste("blacBGalpha",dat,".txt", sep = "")
#filetheta5 <- paste("blacBGbeta",dat,".txt", sep = "")
filetheta6 <- paste("blacBGsigma",dat,".txt", sep = "")
filetheta7 <- paste("blacBGrho",dat,".txt", sep = "")
filetheta8 <- paste("blacBGBsig",dat,".txt", sep = "")

filetheta9 <- paste("KappaState",dat,".txt", sep = "")
filetheta10 <- paste("KappaLambda",dat,".txt", sep = "")

###initial H for missing or no missing value  for  kalman filter
Id = diag(nai)
ha = matrix(c(0),nrow=nai)
H = matrix(c(ha,ha,Id),nrow=nai)
oH = H

H_t <- list()                                         #
for(i in 1:alen){                                     #
  NAs <- which(is.na(ax[i, ]))                        #   
  if (length(NAs) == 0){H_t[[i]] = H}                 #            
  else{                                               #
    H_t[[i]] = H[-NAs, , drop=FALSE]                  #
  }                                                   #
}                                                     #


lik = 0
olik = -Inf #-10^2


#######################################################
# H_t <- list()                                         #
# #
# for(i in 1:alen){                                     #
#   #
#   N <- which(is.na(raw[i, ]))/2                       #   
#   #
#   if (length(N) == 0){H_t[[i]] = H}                   #            
#   #
#   else{                                               #
#     NAs <- N[seq(1, length(N), 2)]                    #
#     H_t[[i]] = H[-NAs, , drop=FALSE]                  #
#   }                                                   #
# }                                                     #


#######################################################







