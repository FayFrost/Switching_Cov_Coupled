## input
## lambda_sw : The switching rates, ou-bm or bm-ou.
## pro_pon_ind: The indicator vector of states, 'pn' potential switch, 'OB' - ou to bm, 'BO' - bm to ou.
## pro_pon_time: The vector of switching times for all switches and sampling points

## output
## the loglikelihood of the states given the lambdas.
lam_lik_cov<-function(lamda_sw, pro_pon_ind, pro_pon_time)
{
  
  all_time = pro_pon_time[ which( pro_pon_ind != 'sp') ]
  all_pon = pro_pon_ind[ which( pro_pon_ind != 'sp') ]  
  pn_inx = which(all_pon=='pn')    # Potential
  OB_inx = which(all_pon=='OB')    # OU to BM 
  BO_inx = which(all_pon=='BO')    # BM to OU

  OB_time <- as.vector(all_time[OB_inx])
  BO_time <- as.vector(all_time[BO_inx])
  pn_time <- as.vector(all_time[pn_inx])
  
  num_pon = all_pon
  num_pon[pn_inx] = 0        # Label 0 for 'pn', -1 for 'OB' and 1 for 'BM'
  num_pon[OB_inx] = -1
  num_pon[BO_inx] = 1 
  num_pon = as.numeric(num_pon)
  
  cum_pon =cumsum(num_pon) + nai
  
  OB_prob <- lookup_lamda(lamda_sw, OB_time)[,1]  /Kappa #sapply(OB_time, lookup_lamda, lamda= lamda_sw)
  BO_prob <- lookup_lamda(lamda_sw, BO_time)[,2 ] /Kappa  
  pn_prob <- lookup_lamda(lamda_sw, pn_time) /Kappa
  
  num_pon[OB_inx] = OB_prob  
  num_pon[BO_inx] = BO_prob  
  num_pon[pn_inx] = 1 - (pn_prob[,2]*(nai-cum_pon[pn_inx])  + pn_prob[,1]*(cum_pon[pn_inx]))  ## needs checking ?
  
  sum(log(num_pon))  
}

